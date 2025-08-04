// signifcant portions of the structure of this file are from a old iteration of the HNSWlib project, a side effect of LIGS being an offshoot of what was intially an attempt to speed up graph construction of HNSWlib

#pragma once
#pragma GCC diagnostic ignored "-fpermissive"


#include <shared_mutex>
#include "visited_list_pool.h"
#include "hnswlib.h"
#include <atomic>
#include <random>
#include <stdlib.h>
#include <assert.h>
#include <unordered_set>
#include <forward_list>
#include <map>
#include <list>
#include <set>
#include <omp.h>
#include <time.h>
#include <bitset>
#include <functional>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <thread>
#include <sys/time.h>
#include <sys/stat.h>
#include <cstdint>
#include <limits>
#include <climits>

#define DISCONNECTED UINT_MAX
#define START_OF_BIN UINT_MAX - 1
#define END_OF_BIN UINT_MAX - 2
#define NO_COND UINT_MAX - 3

//
namespace hnswlib
{
    typedef unsigned int tableint;
    typedef unsigned int linklistsizeint;

    template <typename dist_t>
    class PQbinDoubleGraph : public AlgorithmInterface<dist_t>
    {
    public:
        static const tableint max_update_element_locks = 65536;
        PQbinDoubleGraph(SpaceInterface<dist_t> *s)
        {
        }

        PQbinDoubleGraph(SpaceInterface<dist_t> *s, SpaceInterface<dist_t> *s_sub, const std::string &location, bool nmslib = false, size_t max_elements = 0)
        {
            loadIndex(location, s, s_sub, max_elements);
        }

        PQbinDoubleGraph(SpaceInterface<dist_t> *s, SpaceInterface<dist_t> *s_sub, size_t max_elements, size_t vecdim, size_t divisions, size_t centroids, size_t hyperplane_count, size_t pq_table_no = 0, size_t hyper_table_no = 0, size_t M = 16, size_t ef = 32) : delete_time_indiv(max_elements), insert_time(max_elements)
        {
            max_elements_ = max_elements;

            num_deleted_ = 0;
            num_set_ = 0;
            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();
            data_size_sub_ = s_sub->get_data_size();
            fstdistfunc_sub_ = s_sub->get_dist_func();
            dist_func_param_sub_ = s_sub->get_dist_func_param();
            std::cout << "data_size_sub_: " << data_size_sub_ << "\n";
            std::cout << "data_size_: " << data_size_ << "\n";
            build_time = 0.0;
            delete_time = 0.0;
            total_time = 0.0;
            vecdim = vecdim;

            source_file = 0;
            delete_method = -1;
            update_method = -1;
            delete_selection = -1;
            cutoff_saved = 0.0;
            delete_k = -1;
            bins_initialized = 0;
            delete_chance = -1;
            mode = 0;
            insert_steps = -1;
            iteration = -1;
            M_ = M;
            maxM_ = M_;
            maxM0_ = M_ * 2;

            // file_name = "";
            // gt_address = "";
            // dl_address =  "";
            version = 3;
            point_selection = -1;
            sigma = 0.0;
            bin_calculated = 0;
            bin_standard_deviation = 0.0;
            bin_average = 0.0;
            subdivision = divisions;
            centroid_no = centroids;
            pq_table_number = pq_table_no;
            label_offset_ = data_size_;
            offsetLevel0_ = 0;
            hyperplane_no = hyperplane_count;
            hyperplane_table_number = hyper_table_no;
            total_table_number = hyperplane_table_number + pq_table_number;
            size_links_level0_ = 2 * sizeof(tableint) * total_table_number;
            link_offset_ = data_size_ + sizeof(labeltype) + 2 + sizeof(uint64_t) * total_table_number;
            size_graph_links_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
            size_data_per_element_ = link_offset_ + size_links_level0_ + size_graph_links_; // data -> label -> off bit -> table bin no's

            std::vector<std::vector<std::vector<std::vector<dist_t>>>> vec(pq_table_number, std::vector<std::vector<std::vector<dist_t>>>(subdivision, std::vector<std::vector<dist_t>>(centroid_no, std::vector<dist_t>(centroid_no, 0))));

            pq_dist_storage = vec;

            data_level0_memory_ = (char *)malloc(max_elements_ * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            cur_element_count = 0;

            visited_list_pool_ = new VisitedListPool(1, max_elements);
            std::cout << "pq_table_number: " << pq_table_number << "\n";
            bin_write_mutexes_ = std::vector<std::shared_mutex>(total_table_number);
            bin_lookup_mutexes_ = std::vector<std::mutex>(max_elements);
            // initializations for special treatment of the first node
            std::cout << "centroid_no: " << centroid_no << "\n";
            centroidLists_pq_ = (char *)malloc(data_size_sub_ * subdivision * centroid_no * pq_table_number);
            if (pq_table_number != 0 && centroidLists_pq_ == nullptr)
                throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate bnlists");
            std::cout << "vecdim  *max_elem: " << vecdim * max_elements_ << "\n";
            std::cout << "malloc  centroidLists_pq_: " << data_size_sub_ * subdivision * centroid_no * pq_table_number << "\n";

            std::cout << "hyperplane_no: " << hyperplane_no << "\n";
            centroidLists_hyperplane_ = (char *)malloc(data_size_ * hyperplane_no * hyperplane_table_number);
            if (hyperplane_table_number != 0 && centroidLists_hyperplane_ == nullptr)
                throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate bnlists");
            std::cout << "vecdim  *max_elem: " << vecdim * max_elements_ << "\n";
            std::cout << "malloc  centroidLists_hyperplane_: " << data_size_ * hyperplane_no * hyperplane_table_number << "\n";

            bin_lookup_ = std::vector<std::unordered_map<uint64_t, tableint>>(total_table_number);
            std::cout << "total_table_number: " << total_table_number << "\n";
            for (int i = 0; i < total_table_number; i++)
            {
                // std::unordered_map<uint64_t, std::forward_list<tableint>> *temp = new std::unordered_map<uint64_t, std::forward_list<tableint>>();
                // bin_lookup_.push_back(*temp);
                bin_lookup_[i].reserve(max_elements);
                std::cout << "size of bin lookup: " << bin_lookup_[i].size() << "\n";
            }

            std::cout << "size of bin lookup 0 : " << bin_lookup_[0].size() << "\n";
        }

        size_t max_elements_;
        size_t cur_element_count;
        size_t size_data_per_element_;
        size_t bin_data_per_element_;
        size_t size_graph_links_;
        size_t num_deleted_;
        size_t num_set_;

        // done
        float build_time;
        float delete_time;
        float total_time;
        int subdivision;
        int centroid_no;
        int hyperplane_no;
        int hyperplane_table_number;
        int total_table_number;
        int pq_table_number;
        int vecdim;
        float cutoff_saved;
        // additional
        int source_file;
        int delete_method;
        int update_method;
        int delete_selection;
        int delete_k;
        float delete_chance;
        int mode;
        int insert_steps;
        int iteration;
        // std::string file_name;
        // std::string gt_address;
        // std::string dl_address;
        int version;
        int point_selection;
        float sigma;
        int bin_calculated;
        int bins_initialized;
        float bin_standard_deviation;
        float bin_average;
        int hyperplane;

        VisitedListPool *visited_list_pool_;
        std::mutex cur_element_count_guard_;
        std::vector<std::shared_mutex> bin_write_mutexes_;
        std::vector<std::mutex> bin_lookup_mutexes_;
        size_t offsetData_, offsetLevel0_;

        char *data_level0_memory_;
        char *centroidLists_pq_;
        char *centroidLists_hyperplane_;
        std::vector<float> delete_time_indiv;
        std::vector<float> insert_time;
        std::vector<std::vector<std::vector<std::vector<dist_t>>>> pq_dist_storage;
        size_t label_offset_;
        size_t data_size_;
        size_t link_offset_;
        DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_;
        size_t data_size_sub_;
        DISTFUNC<dist_t> fstdistfunc_sub_;
        void *dist_func_param_sub_;
        std::unordered_map<labeltype, tableint> label_lookup_;
        std::vector<std::unordered_map<uint64_t, tableint>> bin_lookup_;

        size_t M_;
        size_t maxM_;
        size_t maxM0_;
        size_t size_links_level0_;

        struct CompareByFirst
        {
            constexpr bool operator()(std::pair<dist_t, tableint> const &a,
                                      std::pair<dist_t, tableint> const &b) const noexcept
            {
                return a.first < b.first;
            }
        };

        float *fvecs_read(const char *fname, size_t *d_out, size_t *n_out)
        {
            FILE *f = fopen(fname, "r");
            if (!f)
            {
                fprintf(stderr, "could not open %s\n", fname);
                perror("");
                abort();
            }
            int d;
            fread(&d, 1, sizeof(int), f);
            assert((d > 0 && d < 1000000) || !"unreasonable dimension");
            fseek(f, 0, SEEK_SET);
            struct stat st;
            fstat(fileno(f), &st);
            size_t sz = st.st_size;
            assert(sz % ((d + 1) * 4) == 0 || !"weird file size");
            size_t n = sz / ((d + 1) * 4);

            *d_out = d;
            *n_out = n;
            float *x = new float[n * (d + 1)];
            size_t nr = fread(x, sizeof(float), n * (d + 1), f);
            assert(nr == n * (d + 1) || !"could not read whole file");

            // shift array to remove row headers
            for (size_t i = 0; i < n; i++)
                memmove(x + i * d, x + 1 + i * (d + 1), d * sizeof(*x));

            fclose(f);
            // std::cout << "d: " << *d_out << " n: " << *n_out << " x: " << *x << "\n";
            return x;
        }

        double elapsed()
        {
            struct timeval tv;
            gettimeofday(&tv, nullptr);
            return tv.tv_sec + tv.tv_usec * 1e-6;
        }

        ~PQbinDoubleGraph()
        {

            free(data_level0_memory_);
            free(centroidLists_pq_);
            free(centroidLists_hyperplane_);

            delete visited_list_pool_;
        }

        struct Xorshift128Plus
        {
            uint64_t s[2];

            Xorshift128Plus(uint64_t seed)
            {
                s[0] = seed;
                s[1] = seed ^ (uint64_t(4101842887655102017)); // A fixed constant to help with randomness
            }

            uint64_t next()
            {
                uint64_t s1 = s[0];
                const uint64_t s0 = s[1];
                const uint64_t result = s0 + s1;
                s[0] = s0;
                s1 ^= s1 << 23;                          // a
                s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
                return result;
            }
        };

        inline labeltype getExternalLabel(tableint internal_id) const
        {
            labeltype return_label;
            memcpy(&return_label, (data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), sizeof(labeltype));
            return return_label;
        }

        inline void setExternalLabel(tableint internal_id, labeltype label) const
        {
            memcpy((data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), &label, sizeof(labeltype));
        }

        inline labeltype *getExternalLabeLp(tableint internal_id) const
        {
            return (labeltype *)(data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_);
        }

        inline char *getDataByInternalId(tableint internal_id) const
        {
            return (data_level0_memory_ + internal_id * size_data_per_element_);
        }

        linklistsizeint *get_graph_list(tableint internal_id) const
        {
            return (linklistsizeint *)(data_level0_memory_ + internal_id * size_data_per_element_ + link_offset_ + size_links_level0_);
        };

        inline tableint *get_link_prev(tableint internal_id, int table_no) const
        {
            return (tableint *)(data_level0_memory_ + internal_id * size_data_per_element_ + link_offset_ + 2 * sizeof(tableint) * table_no);
        };

        inline tableint *get_link_next(tableint internal_id, int table_no) const
        {
            return (tableint *)(data_level0_memory_ + internal_id * size_data_per_element_ + link_offset_ + 2 * sizeof(tableint) * table_no + sizeof(tableint));
        };

        inline unsigned char *get_binlist(tableint internal_id) const
        {
            //
            return (unsigned char *)(data_level0_memory_ + size_data_per_element_ * internal_id + data_size_ + sizeof(labeltype));
            //
        };

        unsigned short int getListCount(linklistsizeint *ptr) const
        {
            return *((unsigned short int *)ptr);
        }

        inline unsigned char *get_binlist(tableint internal_id, int pq_table_no) const
        {
            //
            return (unsigned char *)(data_level0_memory_ + size_data_per_element_ * internal_id + data_size_ + sizeof(labeltype) + 2 + (pq_table_no * sizeof(u_int64_t)));
            //
        };

        inline char *get_centroidlist(int table_no, int cent_no, int subdiv = 0)
        {
            //
            if (table_no < pq_table_number)
            { // std::cout<< "  table_no:  " << table_no  << " cent_no:  " <<cent_no  << " subdiv: " << subdiv << " " << data_size_sub_ * (subdivision * centroid_no * table_no + centroid_no * subdiv + cent_no) <<"\n";

                return (centroidLists_pq_) + (data_size_sub_ * (subdivision * centroid_no * table_no + centroid_no * subdiv + cent_no));
            }
            else
            { // std::cout<< "  table_no:  " << table_no  << " cent_no:  " <<cent_no  << " hyperplane_no: " << hyperplane_no << " " << data_size_ * (hyperplane_no * table_no + cent_no) <<"\n";

                return (centroidLists_hyperplane_) + (data_size_ * (hyperplane_no * (table_no - pq_table_number) + cent_no));
            }
            //
        };

        inline void add_bin_lookup(tableint internalId, std::vector<std::vector<uint64_t>> bin_val)
        {
            for (int i = 0; i < bin_val.size(); i++)
            {
                add_bin_lookup(internalId, bin_val[i], i);
            }
        };

        inline void add_bin_lookup(tableint internalId, std::vector<uint64_t> bin_val)
        {
            for (int i = 0; i < bin_val.size(); i++)
            {
                add_bin_lookup(internalId, bin_val[i], i);
            }
        };

        inline void add_bin_lookup(tableint internalId, std::vector<uint64_t> bin_val, int pq_table_no)
        {
            for (int i = 0; i < bin_val.size(); i++)
            {
                // std::cout << "at tsi " << pq_table_no << " \n";
                add_bin_lookup(internalId, bin_val[i], pq_table_no);
            }
        };

        inline void add_bin_lookup(tableint internalId, uint64_t bin_val, int pq_table_no) {
            
            // Try to insert or update the value.

                std::unique_lock<std::shared_mutex> insert_lock(bin_write_mutexes_[pq_table_no]);
            auto [it, inserted] = bin_lookup_[pq_table_no].insert({bin_val, internalId});
            insert_lock.unlock();
                      
            if (inserted) {
                //std::cout << "inserted\n";
                //std::cout << "pq_table_no = " << pq_table_no << "\n";
                *(get_link_prev(internalId, pq_table_no)) = START_OF_BIN;
                *(get_link_next(internalId, pq_table_no)) = END_OF_BIN;
            } else {
                    tableint temp = bin_lookup_[pq_table_no][bin_val];

                //std::cout << "start_point = " << temp << "\n";
                //std::cout << "pq_table_no = " << pq_table_no << "\n";
                *(get_link_prev(internalId, pq_table_no)) = START_OF_BIN;
                *(get_link_prev(temp, pq_table_no)) = internalId;
                *(get_link_next(internalId, pq_table_no)) = temp;

                bin_lookup_[pq_table_no][bin_val] = internalId;  // Potential overwrite
            }
        }

       /* inline void add_bin_lookup(tableint internalId, uint64_t bin_val, int pq_table_no) {
            
            // Try to insert or update the value.
            unique_lock<std::shared_mutex> pq_lock(bin_write_mutexes_[pq_table_no]);
            unique_lock<std::shared_mutex> int_lock(bin_lookup_mutexes_[internalId]);
            auto [it, inserted] = bin_lookup_[pq_table_no].insert({bin_val, internalId});
            
            if (inserted) {
                //std::cout << "inserted\n";
                //std::cout << "pq_table_no = " << pq_table_no << "\n";
                *(get_link_prev(internalId, pq_table_no)) = START_OF_BIN;
                *(get_link_next(internalId, pq_table_no)) = END_OF_BIN;
            } else {
                    tableint temp = bin_lookup_[pq_table_no][bin_val];
                unique_lock<std::shared_mutex> temp_lock(bin_lookup_mutexes_[temp]);

                //std::cout << "start_point = " << temp << "\n";
                //std::cout << "pq_table_no = " << pq_table_no << "\n";
                *(get_link_prev(internalId, pq_table_no)) = START_OF_BIN;
                *(get_link_prev(temp, pq_table_no)) = internalId;
                *(get_link_next(internalId, pq_table_no)) = temp;

                bin_lookup_[pq_table_no][bin_val] = internalId;  // Potential overwrite
                temp_lock.unlock();
            }
            int_lock.unlock();
            pq_lock.unlock();
        }*/



        tableint deletePoint(labeltype label)
        {
            // std::cout << "deletePoint\n";
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                // throw std::runtime_error("Label not found markDelete");
                std::cout << "Label not found: " << label << "\n";
                return 0;
            }
            else
            {
                tableint internalId = search->second;
                //std::cout << "label: " << label <<  "internalid: " << internalId << "\n";
                deletePointInternal(internalId);
                //std::cout << "isMarkedDeletedExt: " << isMarkedDeletedExt(label) << ", isLabel: " << isLabel(label) << "\n";
                return internalId;
            }
            // std::cout << "What 3\n";
        }

        inline void deletePointInternal(tableint internalId)
        {
            // std::cout << "deletePointInternal\n";
            for (int i = 0; i < total_table_number; i++)
            {

                delete_bin_lookup(internalId, i);
            }
            markDeletedInternal(internalId);
        }

     /*  inline void delete_bin_lookup(tableint internalId, int i)
        {

                //std::cout << "i: " << i << "\n";

                std::unique_lock<std::shared_mutex> int_lock(bin_lookup_mutexes_[internalId]);
                tableint prev = *get_link_prev(internalId, i);
                tableint next = *get_link_next(internalId, i);
                bool valid_next = (next != END_OF_BIN && next != DISCONNECTED);
                bool start = (prev == START_OF_BIN);
                bool valid_prev = (prev != START_OF_BIN && prev != DISCONNECTED );
                if(valid_next && valid_prev){
                    //std::unique_lock<std::shared_mutex> next_lock(bin_lookup_mutexes_[next]);
                    //std::unique_lock<std::shared_mutex> prev_lock(bin_lookup_mutexes_[prev]);
                    *(get_link_prev(next, i)) = prev; 
                    *(get_link_next(prev, i)) = next;
                }else if(valid_next && !valid_prev){
                        //std::unique_lock<std::shared_mutex> next_lock(bin_lookup_mutexes_[next]);
                        uint64_t bin_label = *get_binlist(internalId, i);
                        std::unique_lock<std::shared_mutex> write_lock(bin_write_mutexes_[i]);
                        bin_lookup_[i][bin_label] = next;
                        *(get_link_prev(next, i)) = prev;
                }else if(!valid_next && valid_prev){
                   // std::unique_lock<std::shared_mutex> prev_lock(bin_lookup_mutexes_[prev]);
                        *(get_link_next(prev, i)) = next;
                }else{
                    uint64_t bin_label = *get_binlist(internalId, i);
                    std::unique_lock<std::shared_mutex> erase_lock(bin_write_mutexes_[i]);
                    bin_lookup_[i].erase(bin_label);
                }
                *(get_link_prev(internalId, i)) = DISCONNECTED;
                *(get_link_next(internalId, i)) = DISCONNECTED;
        };*/

    inline void delete_bin_lookup(tableint internalId, int i)
        {
                //std::cout << "i: " << i << "\n";
                tableint prev = *get_link_prev(internalId, i);
                tableint next = *get_link_next(internalId, i);
                //std::cout << "internalId: " << internalId << "\n";
              // std::lock_guard<std::mutex> int_lock(bin_lookup_mutexes_[internalId]);  //locking
                //std::cout << "prev: " << prev << ", internalId: " << internalId << ", next: " << next << "\n";
                if (next != END_OF_BIN && next != DISCONNECTED )
                { // if not last element
                    //std::cout << "next: " << next << "\n";
                   //std::lock_guard<std::mutex> lock_next(bin_lookup_mutexes_[next]);  //locking
                    if (prev == START_OF_BIN  || prev == DISCONNECTED )
                    { // if first element
                        //std::cout << "next: " << next << "\n";
                        std::unique_lock<std::shared_mutex> write_lock(bin_write_mutexes_[i]);
                        uint64_t bin_label = *get_binlist(internalId, i);
                        if(bin_lookup_[i][bin_label] == internalId){
                            bin_lookup_[i][bin_label] = next;}
                        write_lock.unlock();
                        //std::cout << "bin_lookup_[i][bin_label]: " << bin_lookup_[i][bin_label] << "\n";
                    }else{
                        *(get_link_next(prev, i)) = next;
                    }
                    *(get_link_prev(next, i)) = prev;
                    
                   // lock_next.unlock();//locking
                }
                else if (prev != START_OF_BIN && prev != DISCONNECTED )
                { // if last element but not first element
                    //std::cout << "prev: " << prev << "\n";
                  // std::lock_guard<std::mutex> lock_prev(bin_lookup_mutexes_[prev]); //locking
                    *(get_link_next(prev, i)) = next;
                   // lock_prev.unlock();  //locking
                    //std::cout << "*(get_link_next(prev, i))" << *(get_link_next(prev, i));
                }
                else
                { // if first and last element (i.e. only element)
                    std::unique_lock<std::shared_mutex> erase_lock(bin_write_mutexes_[i]);
                    uint64_t bin_label = *get_binlist(internalId, i);
                    bin_lookup_[i].erase(bin_label);
                    erase_lock.unlock();
                    
                }
                *(get_link_prev(internalId, i)) = DISCONNECTED;
                //*(get_link_next(internalId, i)) = DISCONNECTED;
                //int_lock.unlock();
        };

        /*inline void delete_bin_lookup(tableint internalId, uint64_t bin_val, int pq_table_no)
        {


            std::lock_guard<std::mutex> lock(bin_lookup_mutexes_[pq_table_no]);
            auto it = bin_lookup_[pq_table_no].find(bin_val);
            if (it == bin_lookup_[pq_table_no].end())
            {
                //throw std::runtime_error("bin does not exist");
            }
            else
            {
                it->second.remove(internalId);
            }

            if (it->second.empty())
            {
                bin_lookup_[pq_table_no].erase(it);
            }
        };*/

        /* inline void delete_bin_lookup(tableint internalId, uint64_t bin_val, int pq_table_no)
         {

            // std::cout << "single bin dbl\n";
             // std::lock_guard<std::mutex> lock(bin_lookup_mutexes_[pq_table_no]);
             auto it = bin_lookup_[pq_table_no].find(bin_val);
             if (it == bin_lookup_[pq_table_no].end())
             {
                 // throw std::runtime_error("bin does not exist");
             }
             else
             {
                 auto point = std::find(it->second.begin(), it->second.end(), bin_val);
                 if (point != it->second.end())
                 {
                     if (point == it->second.begin())
                     {
                         // Special case: remove the first element
                         it->second.pop_front();
                     }
                     else
                     {
                         auto prev_point = it->second.before_begin();
                         std::advance(prev_point, std::distance(it->second.begin(), point) - 1);
                         it->second.erase_after(prev_point);
                     }
                 }
                 if (it->second.empty())
                 {
                     bin_lookup_[pq_table_no].erase(it);
                 }
             }

             //std::cout << "What 1\n";
         };*/

        /*
        void build_binList(tableint cur_c)
        {
            // std::cout << "in binlist\n";
            //std::cout << "vecdim: " << vecdim << "\n";
            for(int i = 0; i < )
            normalize_angular(cur_c, vecdim, norm_vec);
            render_binList(cur_c, norm_vec, vecdim);
            norm_vec->~vector();delete_
        }

        void build_binList(labeltype label)
        {
            // std::cout << "in binlist\n";
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                // throw std::runtime_error("Label not found markDelete");
                std::cout << "Label not found: " << label << "\n";
            }
            else
            {
                tableint internalId = search->second;
                build_binList(internalId);
            }
        }
        */

        void saveIndex(const std::string &location)
        {
            std::ofstream output(location, std::ios::binary | std::ios::trunc);
            std::streampos position;
            writeBinaryPOD(output, build_time);
            writeBinaryPOD(output, delete_time);

            writeBinaryPOD(output, offsetLevel0_);
            writeBinaryPOD(output, max_elements_);
            writeBinaryPOD(output, cur_element_count);

            writeBinaryPOD(output, size_data_per_element_);
            writeBinaryPOD(output, bin_data_per_element_);
            writeBinaryPOD(output, data_size_sub_);
            writeBinaryPOD(output, data_size_);
            writeBinaryPOD(output, label_offset_);
            writeBinaryPOD(output, offsetData_);
            writeBinaryPOD(output, subdivision);
            writeBinaryPOD(output, centroid_no);
            writeBinaryPOD(output, hyperplane_no);
            writeBinaryPOD(output, hyperplane_table_number);
            writeBinaryPOD(output, total_table_number);
            writeBinaryPOD(output, pq_table_number);
            writeBinaryPOD(output, cutoff_saved);
            writeBinaryPOD(output, size_links_level0_);
            writeBinaryPOD(output, link_offset_);

            std::cout << "size_data_per_element_: " << size_data_per_element_ << "\n";
            std::cout << "size_links_level0_: " << size_links_level0_ << "\n";
            std::cout << "link_offset_: " << link_offset_ << "\n";
            // std::cout << getConnections(100)[5] << "\n";
            writeBinaryPOD(output, version);
            writeBinaryPOD(output, source_file);
            writeBinaryPOD(output, delete_method);
            writeBinaryPOD(output, update_method);
            writeBinaryPOD(output, delete_selection);
            writeBinaryPOD(output, delete_k);
            writeBinaryPOD(output, delete_chance);
            writeBinaryPOD(output, mode);
            writeBinaryPOD(output, insert_steps);
            writeBinaryPOD(output, iteration);
            writeBinaryPOD(output, point_selection);
            writeBinaryPOD(output, sigma);
            writeBinaryPOD(output, total_time);
            writeBinaryPOD(output, bins_initialized);
            writeBinaryPOD(output, bin_calculated);
            writeBinaryPOD(output, bin_average);
            writeBinaryPOD(output, bin_standard_deviation);

            output.write(data_level0_memory_, cur_element_count * size_data_per_element_);
            output.write(centroidLists_pq_, data_size_sub_ * subdivision * centroid_no * pq_table_number);
            output.write(centroidLists_hyperplane_, data_size_ * hyperplane_no * hyperplane_table_number);

            for (size_t i = 0; i < cur_element_count; i++)
            {
                // output2 << std::to_string(search_not_deleted[i]) <<  std::endl;
                // output2 << std::to_string(search_deleted[i]) <<  std::endl;
                writeBinaryPOD(output, delete_time_indiv[i]);
                writeBinaryPOD(output, insert_time[i]);
            }

            for (size_t i = 0; i < pq_table_number; i++)
            {
                for (size_t j = 0; j < subdivision; j++)
                {
                    for (size_t k = 0; k < centroid_no; k++)
                    {
                        for (size_t l = 0; l < centroid_no; l++)
                        {

                            writeBinaryPOD(output, pq_dist_storage[i][j][k][l]);
                        }
                    }
                }
            }
            
            // output.write(binLists_, vecdim * cur_element_count);
            output.close();
            // output2.close();

            std::ofstream output2(location + "_add2", std::ios::binary | std::ios::trunc);
            writeBinaryPOD(output2, subdivision);
            writeBinaryPOD(output2, centroid_no);
            writeBinaryPOD(output2, hyperplane_no);
            writeBinaryPOD(output2, hyperplane_table_number);
            writeBinaryPOD(output2, total_table_number);
            writeBinaryPOD(output2, pq_table_number);
            output2.close();
        }

        void loadIndex(const std::string &location, SpaceInterface<dist_t> *s, SpaceInterface<dist_t> *s_sub, size_t max_elements_i = 0)
        {
            std::ifstream input(location, std::ios::binary);

            if (!input.is_open())
                throw std::runtime_error("Cannot open file");

            // get file size:
            input.seekg(0, input.end);
            std::streampos total_filesize = input.tellg();
            input.seekg(0, input.beg);

            // done
            readBinaryPOD(input, build_time);
            readBinaryPOD(input, delete_time);

            readBinaryPOD(input, offsetLevel0_);
            readBinaryPOD(input, max_elements_);
            readBinaryPOD(input, cur_element_count);

            size_t max_elements = max_elements_i;
            if (max_elements < cur_element_count)
                max_elements = max_elements_;
            max_elements_ = max_elements;

            readBinaryPOD(input, size_data_per_element_);
            std::cout << "size_data_per_element_: " << size_data_per_element_ << "\n";
            readBinaryPOD(input, bin_data_per_element_);
            readBinaryPOD(input, data_size_sub_);
            readBinaryPOD(input, data_size_);
            readBinaryPOD(input, label_offset_);
            readBinaryPOD(input, offsetData_);
            readBinaryPOD(input, subdivision);
            readBinaryPOD(input, centroid_no);
            readBinaryPOD(input, hyperplane_no);
            readBinaryPOD(input, hyperplane_table_number);
            readBinaryPOD(input, total_table_number);
            readBinaryPOD(input, pq_table_number);
            readBinaryPOD(input, cutoff_saved);
            readBinaryPOD(input, size_links_level0_);
            readBinaryPOD(input, link_offset_);
            std::cout << "size_links_level0_: " << size_links_level0_ << "\n";
            // std::cout << "test 1\n";
            readBinaryPOD(input, version);
            readBinaryPOD(input, source_file);
            readBinaryPOD(input, delete_method);
            readBinaryPOD(input, update_method);
            readBinaryPOD(input, delete_selection);
            readBinaryPOD(input, delete_k);
            readBinaryPOD(input, delete_chance);
            readBinaryPOD(input, mode);
            readBinaryPOD(input, insert_steps);
            readBinaryPOD(input, iteration);
            readBinaryPOD(input, point_selection);
            readBinaryPOD(input, sigma);
            readBinaryPOD(input, total_time);
            readBinaryPOD(input, bins_initialized);
            readBinaryPOD(input, bin_calculated);
            readBinaryPOD(input, bin_average);
            readBinaryPOD(input, bin_standard_deviation);
            // std::cout << "test 2\n";

            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();
            data_size_sub_ = s_sub->get_data_size();
            fstdistfunc_sub_ = s_sub->get_dist_func();
            dist_func_param_sub_ = s_sub->get_dist_func_param();
            // std::cout << "test 3\n";

            auto pos = input.tellg();

            input.seekg(pos, input.beg);
            // std::cout << "test 4\n";


            bin_write_mutexes_ = std::vector<std::shared_mutex>(total_table_number);
            bin_lookup_mutexes_ = std::vector<std::mutex>(max_elements);
            data_level0_memory_ = (char *)malloc(max_elements * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");
            centroidLists_pq_ = (char *)malloc(data_size_sub_ * subdivision * centroid_no * pq_table_number);
            if (centroidLists_pq_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate centroidLists_ ");
            centroidLists_hyperplane_ = (char *)malloc(data_size_ * hyperplane_no * hyperplane_table_number);
            if (centroidLists_hyperplane_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate centroidLists_ ");
            input.read(data_level0_memory_, cur_element_count * size_data_per_element_);
            input.read(centroidLists_pq_, data_size_sub_ * subdivision * centroid_no * pq_table_number);
            input.read(centroidLists_hyperplane_, data_size_ * hyperplane_no * hyperplane_table_number);

            visited_list_pool_ = new VisitedListPool(1, max_elements);
            // std::cout << "test 6\n";

            for (size_t i = 0; i < cur_element_count; i++)
            {
                label_lookup_[getExternalLabel(i)] = i;
                // speedup unmarkBinSetInternal(i);
            }
            /*
            for (size_t i = 0; i < cur_element_count; i++)
            {
                if (isMarkedDeleted(i))
                {
                    num_deleted_ += 1;
                }
                //if (isMarkedBinSet(i))
                {
                    num_set_ += 1;
                    uint64_t bin_val = getBin(i);
                    add_bin_lookup(i, bin_val);
                }
            }
            */

            delete_time_indiv = std::vector<float>(max_elements);
            insert_time = std::vector<float>(max_elements);
            // std::cout << "test 7\n";

            for (size_t i = 0; i < cur_element_count; i++)
            {

                readBinaryPOD(input, delete_time_indiv[i]);
                readBinaryPOD(input, insert_time[i]);
                /*std::cout << "i: " << i << " \n";
                std::string temp;
                input_add2 >> temp;
                std::cout << "temp: " << std::stoi(temp) << " \n";
                search_not_deleted[i] = std::stoi(temp);
                input_add2 >> search_deleted[i];
                search_deleted[i] = std::stoi(temp);*/
            }

            // std::cout << "test 8\n";
            std::vector<std::vector<std::vector<std::vector<dist_t>>>> vec(pq_table_number, std::vector<std::vector<std::vector<dist_t>>>(subdivision, std::vector<std::vector<dist_t>>(centroid_no, std::vector<dist_t>(centroid_no, 0))));

            pq_dist_storage = vec;
            // std::cout << "test 8.\n";
            for (size_t i = 0; i < pq_table_number; i++)
            {
                for (size_t j = 0; j < subdivision; j++)
                {
                    for (size_t k = 0; k < centroid_no; k++)
                    {
                        for (size_t l = 0; l < centroid_no; l++)
                        {

                            readBinaryPOD(input, pq_dist_storage[i][j][k][l]);
                        }
                    }
                }
            }

            // std::cout << "test 9\n";
            bin_lookup_ = std::vector<std::unordered_map<uint64_t, tableint>>(total_table_number);
            std::cout << "total_table_number: " << total_table_number << "\n";
            // for (int i = 0; i < total_table_number; i++)
            // {
            //     std::unordered_map<uint64_t, std::forward_list<tableint>> *temp = new std::unordered_map<uint64_t, std::forward_list<tableint>>();
            //     bin_lookup_.push_back(*temp);
            //     std::cout << "size of bin lookup: " << bin_lookup_[i].size() << "\n";
            // }
            std::cout << "size of bin lookup 0 : " << bin_lookup_[0].size() << "\n";
            std::cout << "cur_element_count : " << cur_element_count << "\n";

            // std::cout << "test 8\n";

            // speedup initializeBins();

            for (size_t i = 0; i < cur_element_count; i++)
            {
                std::vector<u_int64_t> bin_val(total_table_number);
                unsigned char *bin_cur = (get_binlist(i)) + 2;
                for (int j = 0; j < total_table_number; j++)
                {

                    memcpy(&(bin_val[j]), bin_cur + j * sizeof(u_int64_t), sizeof(u_int64_t));

                    // std::cout << bin_val[j] << "\n";
                }
                add_bin_lookup(i, bin_val);
                // speedup unmarkBinSetInternal(i);
            }

            // std::cout << "test 9\n";

            std::cout << "size of bin lookup 0 : " << bin_lookup_[0].size() << "\n";
            input.close();

            // input_add2.close();

            return;
        }

        bool isLabel(labeltype label) const
        {
            tableint label_c;
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end() || isMarkedDeleted(search->second))
            {
                return false;
            }
            return true;
        }

        template <typename data_t>
        std::vector<data_t> getDataByLabel(labeltype label) const
        {
            tableint label_c;
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end() || isMarkedDeleted(search->second))
            {
                throw std::runtime_error("Label not found getDatabyLabel");
            }
            label_c = search->second;

            char *data_ptrv = getDataByInternalId(label_c);
            size_t dim = *((size_t *)dist_func_param_);
            std::vector<data_t> data;
            data_t *data_ptr = (data_t *)data_ptrv;
            for (int i = 0; i < dim; i++)
            {
                data.push_back(*data_ptr);
                data_ptr += 1;
            }
            return data;
        }

        template <bool ignore_delete = false>
        inline char *getDataByLabelChar(labeltype label) const
        {
            tableint label_c;
            auto search = label_lookup_.find(label);
            if (ignore_delete)
            {
                if (search == label_lookup_.end())
                {
                    throw std::runtime_error("Label not found getDatabyLabel 1");
                }
            }
            else
            {
                if (search == label_lookup_.end() || isMarkedDeleted(search->second))
                {
                    throw std::runtime_error("Label not found getDatabyLabel 2");
                }
            }

            label_c = search->second;

            char *data_ptrv = getDataByInternalId(label_c);
            return data_ptrv;
        }

        static const unsigned char NORMAL_MARK = 0x00;
        static const unsigned char DELETE_MARK = 0x01;
        static const unsigned char UNSET_MARK = 0x00;
        static const unsigned char SET_MARK = 0x01;
        // static const unsigned char REUSE_MARK = 0x10;
        /**
         * Marks an element with the given label deleted, does NOT really change the current graph.
         * @param label
         */
        tableint markDelete(labeltype label)
        {
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                // throw std::runtime_error("Label not found markDelete");
                std::cout << "Label not found: " << label << "\n";
                return 0;
            }
            else
            {
                tableint internalId = search->second;
                markDeletedInternal(internalId);
                return internalId;
            }
        }

        /**
         * Uses the first 8 bits of the memory for the linked list to store the mark,
         * whereas maxM0_ has to be limited to the lower 24 bits, however, still large enough in almost all cases.
         * @param internalId
         */

        bool isMarkedBinSet(tableint internalId)
        {
            //
            unsigned char *ll_cur = (get_binlist(internalId)) + 1;
            return *ll_cur & SET_MARK;
            //
        };

        bool markBinSetInternal(tableint internalId)
        {
            // std::cout << "cur_element_count " << cur_element_count << "\n";
            // std::cout << "internalId " << internalId << "\n";

            //
            if (isMarkedBinSet(internalId))
            {
                //   std::cout << "is marked bin set "
                //           << "\n";

                return false;
            }
            else
            {
                //  std::cout << "is not marked bin set "
                //         << "\n";
                unsigned char *ll_cur = (get_binlist(internalId)) + 1;
                // std::cout << (uint8_t)ll_cur[0] <<"\n";
                *ll_cur |= SET_MARK;

                // std::cout << (uint8_t)ll_cur[0] <<"\n";
                num_set_++;

                return true;
            }
            //
        };

        void unmarkBinSetInternal(tableint internalId)
        {

            // assert(internalId < cur_element_count);
            if (isMarkedBinSet(internalId))
            {
                unsigned char *ll_cur = (get_binlist(internalId)) + 1;
                *ll_cur &= ~SET_MARK;
                num_set_ -= 1;
            }
            else
            {
                throw std::runtime_error("The requested to unset element is not set");
            }
        }

        int myPow(uint64_t x, uint64_t p)
        {
            if (p == 0)
                return 1;
            if (p == 1)
                return x;

            uint64_t tmp = myPow(x, p / 2);
            if (p % 2 == 0)
                return tmp * tmp;
            else
                return x * tmp * tmp;
        }

        void unmarkBinSet(labeltype label)
        {
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                throw std::runtime_error("Label not found unmarkBinSet");
            }
            tableint internalId = search->second;
            unmarkBinSetInternal(internalId);
        }

        void markBinSet(labeltype label)
        {
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                throw std::runtime_error("Label not found markBinSet");
            }
            tableint internalId = search->second;
            markBinSetInternal(internalId);
        }

        void markDeletedInternal(tableint internalId)
        {
            assert(internalId < cur_element_count);
            if (!isMarkedDeleted(internalId))
            {
                unsigned char *ll_cur = (get_binlist(internalId));
                *ll_cur |= DELETE_MARK;
                num_deleted_ += 1;
            }
            else
            {
                std::cout << "The requested to delete element is already deleted: " << getExternalLabel(internalId) << "\n";
                // throw std::runtime_error("The requested to delete element is already deleted");
            }
        }

        /**
         * Remove the deleted mark of the node, does NOT really change the current graph.
         * @param label
         */
        void unmarkDelete(labeltype label)
        {
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                throw std::runtime_error("Label not found unmarkDelete");
            }
            tableint internalId = search->second;
            unmarkDeletedInternal(internalId);
        }

        void markNormalInternal(tableint internalId)
        {

            assert(internalId < cur_element_count);
            if (isMarkedDeleted(internalId))
            {
                unsigned char *ll_cur = (get_binlist(internalId));
                *ll_cur = NORMAL_MARK;
            }
        }

        /**
         * Remove the deleted mark of the node.
         * @param internalId
         */
        void unmarkDeletedInternal(tableint internalId)
        {

            assert(internalId < cur_element_count);
            if (isMarkedDeleted(internalId))
            {
                unsigned char *ll_cur = (get_binlist(internalId));
                ;
                *ll_cur &= ~DELETE_MARK;
                num_deleted_ -= 1;
            }
            else
            {
                throw std::runtime_error("The requested to undelete element is not deleted");
            }
        }

        void normalize_angular(labeltype label, int vecdim)
        {

            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end())
            {
                throw std::runtime_error("Label not found normalize_angular");
            }
            tableint internalId = search->second;
            normalize_angular(internalId, vecdim);
        }

        void normalize_angular(tableint internalId, int vecdim)
        {
            char *mass = getDataByInternalId(internalId);
            std::vector<float> norm_vec(vecdim);
            float total = 0;
            for (int i = 0; i < vecdim * 4; i = i + 4)
            {
                char temp[4];
                for (int j = 0; j < 4; j++)
                {
                    temp[j] = mass[i + j];
                }
                float *temp_float = (float *)temp;
                // std::cout << "temp_float at " <<  i/4 << ": " << *temp_float << "\n";
                norm_vec[i / 4] = *temp_float;
                total = total + norm_vec[i / 4] * norm_vec[i / 4];
                // std::cout << "total: " << total << "\n";
            }
            float total_sqrt = sqrt(total);
            for (int i = 0; i < vecdim; i++)
            {
                char temp[4];
                // std::cout << "norm_vec[" << i << "]: " <<  norm_vec[i] << "\n";
                norm_vec[i] = norm_vec[i] / total_sqrt;

                // std::cout << "norm_vec[" << i << "]: " <<  norm_vec[i] << "\n";
                memcpy(temp, &(norm_vec[i]), sizeof(float));
                for (int j = 0; j < 4; j++)
                {
                    mass[i * 4 + j] = temp[j];
                }
                /*
                float *temp_float = (float *)temp;
                std::cout << "temp_float at " <<  i << ": " << *temp_float << "\n";
                char *mass2 = getDataByInternalId(internalId);
                char temp2[4];
                for (int j = 0; j < 4; j++)
                {
                    temp2[j] = mass2[i * 4 + j];
                }
                 float *temp_float2 = (float *)temp2;
                std::cout << "temp_float2 at " <<  i << ": " << *temp_float2 << "\n";*/
            }
        }
        /**
         * Checks the first 8 bits of the memory to see if the element is marked deleted.
         * @param internalId
         * @return
         */
        bool isMarkedDeleted(tableint internalId) const
        {
            unsigned char *ll_cur = (get_binlist(internalId));
            return *ll_cur & DELETE_MARK;
        }

        bool isMarkedDeletedExt(labeltype externalID) const
        {
            auto search = label_lookup_.find(externalID);
            if (search == label_lookup_.end())
            {
                throw std::runtime_error("Label not found isMarkedDeletedExt");
            }
            return isMarkedDeleted(search->second);
        }

        void updatePointLabel(const void *dataPoint, labeltype label, bool cutoff_applied = false, float cutoff = 1.0)
        {
            // TODO : add check for lable existing
            auto search = label_lookup_.find(label);
            if (search != label_lookup_.end())
            {
                tableint existingInternalId = search->second;
                updatePoint<true>(dataPoint, existingInternalId);
            }
            else
            {
                std::cout << "Point does not exist - using addPoint";
                addPoint<true>(dataPoint, label, cutoff_applied, cutoff);
            }
        }

        tableint updatePoint(tableint internalId)
        {
            tableint temp = updatePoint(getDataByInternalId(internalId), internalId);
            return temp;
        }

        std::vector<uint64_t> getBin(tableint internalId)
        {
            std::vector<uint64_t> temp_v;
            for (int i = 0; i < total_table_number; i++)
            {
                int temp = getBin(internalId, i);
                temp_v.push_back(temp);
            }
            return temp_v;
        }

        uint64_t getBin(tableint internalId, int pq_table_no)
        {
            unsigned char *bin_cur = (get_binlist(internalId, pq_table_no));
            char temp[8];
            for (int j = 0; j < 8; j++)
            {
                temp[j] = bin_cur[j];
            }
            uint64_t *temp_ui64 = (uint64_t *)temp;
            return *temp_ui64;
        }

        uint64_t calculateBin(const void *dataPoint, int pq_bin_no)
        {
            // std::cout << "subdiv: " << subdivision <<"\n";
            // std::cout << "centroid_no: " << centroid_no <<"\n";
            uint64_t bin_val = 0;
            if (pq_bin_no < pq_table_number)
            {
                for (int i = 0; i < subdivision; i++)
                {
                    bin_val = bin_val * centroid_no;
                    char *temp = ((char *)dataPoint) + i * data_size_sub_;
                    float dist = std::numeric_limits<float>::max();

                    int chosen_centroid;
                    for (int j = 0; j < centroid_no; j++)
                    {
                        char *bin_cur = get_centroidlist(pq_bin_no, j, i);
                        float cent_dist = fstdistfunc_sub_(temp, bin_cur, dist_func_param_sub_);
                        if (cent_dist < dist)
                        {
                            dist = cent_dist;
                            chosen_centroid = j;
                        }
                    }
                    bin_val = bin_val + chosen_centroid;
                    // std::cout << "subdivision: " << i << ", centroid_no: " << chosen_centroid << ", dist: " << dist << "\n";
                }

                return bin_val;
            }
            else
            {
                uint64_t bin_val = 0;
                for (int j = 0; j < hyperplane_no; j++)
                {
                    float *bin_cur = (float *)get_centroidlist(pq_bin_no, j, 0);
                    float sum_ = 0;
                    for (int i = 0; i < vecdim; i++)
                    {
                        sum_ = bin_cur[i] * ((float *)dataPoint)[i] + sum_;

                        // std::cout << "j: " << j << ", i: " << i << " bin_cur: " << bin_cur[i] << ", dataPoint: " << ((float *)dataPoint)[i] << ", sum_: " << sum_ << "\n";
                    }
                    bin_val = bin_val * 2;
                    if (sum_ > 0)
                    {
                        bin_val++;
                    }
                }
                // std::cout << "subdivision: " << i << ", centroid_no: " << chosen_centroid << ", dist: " << dist << "\n";
                // std::cout << "bin_val: " << bin_val << "\n";
                return bin_val;
            }
        }

        std::vector<uint64_t> calculateBin(const void *dataPoint)
        {
            std::vector<uint64_t> bins_all_sets(total_table_number);
            {
                for (int i = 0; i < total_table_number; i++)
                {
                    bins_all_sets[i] = calculateBin(dataPoint, i);
                }
            }

            return bins_all_sets;
        }

        uint64_t assignBin(tableint internalId)
        {

            uint64_t newbin = (calculateBin(getDataByInternalId(internalId)))[0];

            setBin(internalId, newbin);
            if (!isMarkedBinSet(internalId))
            {
                markBinSetInternal(internalId); // done
            };
            return newbin;
        }

        uint64_t assignBinLabel(labeltype label)
        {

            auto search = label_lookup_.find(label);
            if (search != label_lookup_.end())
            {
                tableint existingInternalId = search->second;
                return assignBin(existingInternalId);
            }
            else
            {
                throw std::runtime_error("Assign bin failed - label not found");
            }
        }

        std::priority_queue<std::pair<dist_t, labeltype>> searchKnn(const void *query_data, size_t k)
        {
            int temp = 0.0;
            int temp2 = 0.0;
            int temp3 = 0.0;
            return searchKnn_Hamming(query_data, k, 1, temp, temp2, temp3, nullptr, false);
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayerST(tableint ep_id, const void *data_point, size_t ef, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, int cutoff_type = 0, float cutoff = 2.0)
        { // done
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;
    //std::cout <<" base layer st\n";
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

            dist_t lowerBound;
            dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            lowerBound = dist;
            top_candidates.emplace(dist, ep_id);
            candidate_set.emplace(-dist, ep_id);

            visited_array[ep_id] = visited_array_tag;
            // std::cout << "before candiate_set\n";
            std::set<u_int64_t> checked_bins;

            while (!candidate_set.empty())
            {

                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();

                if ((-current_node_pair.first) > lowerBound && (top_candidates.size() == ef))
                {
                    break;
                }
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;

                u_int64_t *data = (u_int64_t *)(get_binlist(current_node_id) + 2);
                //                bool cur_node_deleted = isMarkedDeleted(current_node_id);

                // #ifdef USE_SSE
                //                 _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
                //                 _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
                //                 _mm_prefetch(data_level0_memory_ + (*(data + 1)) * size_data_per_element_ + offsetData_, _MM_HINT_T0);
                //                 _mm_prefetch((char *)(data + 2), _MM_HINT_T0);
                // #endif

                for (int pqtn = 0; pqtn < total_table_number; pqtn++)
                {

                    // std::cout << "pqtn" << pqtn  << " \n";
                    u_int64_t *candidate_bin = (u_int64_t *)get_binlist(current_node_id, pqtn);

                    checked_bins.insert(*candidate_bin * 100 + pqtn);

                    int point_tracker = 0;
                    // std::cout << "candid.bin: " << *candidate_bin << "\n";
                    if (bin_lookup_[pqtn].find(*candidate_bin) != bin_lookup_[pqtn].end())
                    {
                        
                        tableint cand_ = bin_lookup_[pqtn][*candidate_bin];
                        //std::cout << "cand_: " << cand_ << "\n";
                        // for (tableint &cand_ : bin_lookup_[pqtn][*candidate_bin])
                        while (cand_ < NO_COND)
                        {

                            // (checked_points_no)++;
                            // point_tracker++;

                            // #ifdef USE_SSE
                            //                     _mm_prefetch((char *)(visited_array + *(data + j + 1)), _MM_HINT_T0);
                            //                     _mm_prefetch(data_level0_memory_ + (*(data + j + 1)) * size_data_per_element_ + offsetData_,
                            //                                  _MM_HINT_T0); ////////////
                            // #endif

                            if (!(visited_array[cand_] == visited_array_tag))
                            {

                                visited_array[cand_] = visited_array_tag;
                                //(visited_points_no)++;
                                char *currObj1 = (getDataByInternalId(cand_));

                                dist_t dist = fstdistfunc_(data_point, currObj1, dist_func_param_);

                                if (top_candidates.size() < ef || lowerBound > dist)
                                {

                                    candidate_set.emplace(-dist, cand_);
                                    // #ifdef USE_SSE
                                    //                             _mm_prefetch(data_level0_memory_ + candidate_set.top().second * size_data_per_element_ +
                                    //                                              offsetLevel0_, ///////////
                                    //                                          _MM_HINT_T0);      ////////////////////////
                                    // #endif

                                    if (!isMarkedDeleted(cand_))
                                    {
                                        top_candidates.emplace(dist, cand_);
                                    }

                                    if (top_candidates.size() > ef)
                                    {
                                        top_candidates.pop();
                                    }

                                    if (!top_candidates.empty())
                                    {
                                        lowerBound = top_candidates.top().first;
                                    }
                                }
                            }

                            //std::cout << "cand_: " << cand_;
                            cand_ = *get_link_next(cand_, pqtn);
                            //std::cout << ", next cand_: " << cand_ << "\n";
                        }
                    }
                }
            }

            visited_list_pool_->releaseVisitedList(vl);
            checked_bins_no = checked_bins.size();
            return top_candidates;
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayerST_inf(tableint ep_id, const void *data_point, size_t ef, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, int cutoff_type = 1, float cutoff = 1000.0)
        { // done
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;
            dist_t cutoff_1_bound;
            dist_t cutoff_2_bound;
            dist_t lowerBound;
            cutoff = 1.3;
            dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            lowerBound = dist;
            cutoff_1_bound = dist;
            cutoff_2_bound = dist;
            top_candidates.emplace(dist, ep_id);
            candidate_set.emplace(-dist, ep_id);
            visited_array[ep_id] = visited_array_tag;
            // std::cout << "before candiate_set\n";
            std::set<u_int64_t> checked_bins;

            while (!candidate_set.empty())
            {

                if (cutoff_1_bound > -candidate_set.top().first)
                {
                    cutoff_1_bound = -candidate_set.top().first;
                }
                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();
                if ((-current_node_pair.first) > lowerBound && (top_candidates.size() == ef))
                {
                    break;
                }
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;

                u_int64_t *data = (u_int64_t *)(get_binlist(current_node_id) + 2);
                //                bool cur_node_deleted = isMarkedDeleted(current_node_id);

                // #ifdef USE_SSE
                //                 _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
                //                 _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
                //                 _mm_prefetch(data_level0_memory_ + (*(data + 1)) * size_data_per_element_ + offsetData_, _MM_HINT_T0);
                //                 _mm_prefetch((char *)(data + 2), _MM_HINT_T0);
                // #endif

                for (int pqtn = 0; pqtn < total_table_number; pqtn++)
                {

                    // std::cout << "pqtn" << pqtn  << " \n";
                    u_int64_t *candidate_bin = (u_int64_t *)get_binlist(current_node_id, pqtn);

                    checked_bins.insert(*candidate_bin * 100 + pqtn);

                    int point_tracker = 0;
                    // std::cout << "candid.bin: " << *candidate_bin << "\n";
                    if (bin_lookup_[pqtn].find(*candidate_bin) != bin_lookup_[pqtn].end())
                    {
                        tableint cand_ = START_OF_BIN;
                        cand_ = bin_lookup_[pqtn][*candidate_bin];

                        while (cand_< NO_COND)
                        {
                            //(checked_points_no)++;
                            // point_tracker++;

                            // #ifdef USE_SSE
                            //                     _mm_prefetch((char *)(visited_array + *(data + j + 1)), _MM_HINT_T0);
                            //                     _mm_prefetch(data_level0_memory_ + (*(data + j + 1)) * size_data_per_element_ + offsetData_,
                            //                                  _MM_HINT_T0); ////////////
                            // #endif

                            if (!(visited_array[cand_] == visited_array_tag))
                            {

                                visited_array[cand_] = visited_array_tag;
                                //(visited_points_no)++;
                                char *currObj1 = (getDataByInternalId(cand_));

                                dist_t dist = fstdistfunc_(data_point, currObj1, dist_func_param_);

                                if (top_candidates.size() < ef || lowerBound > dist)
                                {

                                    // std::cout  << "cutoff_:" << cutoff_1_bound << " dist:" << dist << " -dist: " << -dist << "\n";
                                    if (dist < cutoff_1_bound * cutoff)
                                    {
                                        candidate_set.emplace(-dist, cand_);
                                    }
                                    // #ifdef USE_SSE
                                    //                             _mm_prefetch(data_level0_memory_ + candidate_set.top().second * size_data_per_element_ +
                                    //                                              offsetLevel0_, ///////////
                                    //                                          _MM_HINT_T0);      ////////////////////////
                                    // #endif

                                    // note:: the closer one is earlier for candidate_set
                                    if (!isMarkedDeleted(cand_))
                                    {
                                        top_candidates.emplace(dist, cand_);
                                    }

                                    if (top_candidates.size() > ef)
                                    {
                                        top_candidates.pop();
                                    }

                                    if (!top_candidates.empty())
                                    {
                                        lowerBound = top_candidates.top().first;
                                        if (cutoff_1_bound > dist)
                                        {
                                            cutoff_1_bound = dist;
                                        }
                                    }
                                }
                            }
                            cand_ = *get_link_next(cand_, pqtn);
                        }
                        /*#pragma omp critical
                                            {
                                                visited_bins_pop->at(point_tracker)++;
                                            }*/
                    }
                }
            }

            visited_list_pool_->releaseVisitedList(vl);
            checked_bins_no = checked_bins.size();
            return top_candidates;
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayerST_cutoff(tableint ep_id, const void *data_point, size_t ef, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, int cutoff_type, float cutoff)
        { // done
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

            dist_t lowerBound;
            dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            lowerBound = dist;
            top_candidates.emplace(dist, ep_id);
            candidate_set.emplace(-dist, ep_id);

            visited_array[ep_id] = visited_array_tag;
            // std::cout << "before candiate_set\n";
            std::set<u_int64_t> checked_bins;

            while (!candidate_set.empty())
            {

                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();

                if ((-current_node_pair.first) > lowerBound && (top_candidates.size() == ef))
                {
                    break;
                }
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;

                u_int64_t *data = (u_int64_t *)(get_binlist(current_node_id) + 2);
                //                bool cur_node_deleted = isMarkedDeleted(current_node_id);

                // #ifdef USE_SSE
                //                 _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
                //                 _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
                //                 _mm_prefetch(data_level0_memory_ + (*(data + 1)) * size_data_per_element_ + offsetData_, _MM_HINT_T0);
                //                 _mm_prefetch((char *)(data + 2), _MM_HINT_T0);
                // #endif

                int point_tracker = 0;
                // std::cout << "candid.bin: " << *candidate_bin << "\n";
                std::vector<std::pair<dist_t, tableint>> temp_bin_candidates;
                temp_bin_candidates.reserve(ef * total_table_number);

                for (int pqtn = 0; pqtn < total_table_number; pqtn++)
                {

                    // std::cout << "pqtn" << pqtn  << " \n";
                    u_int64_t *candidate_bin = (u_int64_t *)get_binlist(current_node_id, pqtn);

                    checked_bins.insert(*candidate_bin * 100 + pqtn);

                    if (bin_lookup_[pqtn].find(*candidate_bin) != bin_lookup_[pqtn].end())
                    {
                        tableint cand_ = START_OF_BIN;
                        cand_ = bin_lookup_[pqtn][*candidate_bin];

                        // for (tableint &cand_ : bin_lookup_[pqtn][*candidate_bin])
                        while (cand_< NO_COND)
                        {

                            char *currObj1 = (getDataByInternalId(cand_));

                            dist_t dist = fstdistfunc_(data_point, currObj1, dist_func_param_);
                            temp_bin_candidates.push_back(std::make_pair(dist, cand_));

                            cand_ = *get_link_next(cand_, pqtn);
                        }
                    }
                }
                {

                    for (auto cand_pair : temp_bin_candidates)
                    {
                        tableint cand_ = cand_pair.second;
                        tableint dist = cand_pair.first;
                        //(checked_points_no)++;
                        // point_tracker++;

                        // #ifdef USE_SSE
                        //                     _mm_prefetch((char *)(visited_array + *(data + j + 1)), _MM_HINT_T0);
                        //                     _mm_prefetch(data_level0_memory_ + (*(data + j + 1)) * size_data_per_element_ + offsetData_,
                        //                                  _MM_HINT_T0); ////////////
                        // #endif

                        if (!(visited_array[cand_] == visited_array_tag))
                        {

                            visited_array[cand_] = visited_array_tag;
                            //(visited_points_no)++;

                            if (top_candidates.size() < ef || lowerBound > dist)
                            {

                                candidate_set.emplace(-dist, cand_);
                                // #ifdef USE_SSE
                                //                             _mm_prefetch(data_level0_memory_ + candidate_set.top().second * size_data_per_element_ +
                                //                                              offsetLevel0_, ///////////
                                //                                          _MM_HINT_T0);      ////////////////////////
                                // #endif

                                if (!isMarkedDeleted(cand_))
                                {
                                    top_candidates.emplace(dist, cand_);
                                }

                                if (top_candidates.size() > ef)
                                {
                                    top_candidates.pop();
                                }

                                if (!top_candidates.empty())
                                {
                                    lowerBound = top_candidates.top().first;
                                }
                            }
                        }
                    }
                    // #pragma omp critical
                    /*{
                         visited_bins_pop->at(point_tracker)++;
                     }*/
                }
            }

            visited_list_pool_->releaseVisitedList(vl);
            checked_bins_no = checked_bins.size();
            return top_candidates;
        }

        void render_hamming(std::vector<std::pair<int, uint64_t>> *bins_to_check_pair, std::vector<std::set<uint64_t>> *bins_to_check, std::vector<std::set<uint64_t>> *checked_bins)
        {

            std::vector<std::pair<int, uint64_t>> bins_to_check_pair_temp;
            bins_to_check_pair_temp.reserve(total_table_number * pq_table_number * centroid_no);
            for (int m = 0; m < total_table_number; m++)
            {
                std::set<uint64_t> bins_to_check_update;
                if (m < pq_table_number)
                {
                    for (auto itr : (*bins_to_check)[m])
                    {
                        // std::bitset<64> itr_bin(itr);
                        // std::cout << "itr_bin: " << itr_bin << "\n";
                        // std::cout << "btcu: \n";
                        int bins[subdivision];
                        uint64_t temp = itr;
                        for (int i = 0; i < subdivision; i++)
                        {
                            bins[i] = temp % centroid_no;
                            temp = (temp) / centroid_no;
                        }

                        for (int i = 0; i < subdivision; i++)
                        {
                            for (int j = 0; j < centroid_no; j++)
                            {
                                // celcing bins to check
                                uint64_t bin_to_add = 0;

                                for (int k = subdivision - 1; k >= 0; k--)
                                {
                                    bin_to_add = bin_to_add * centroid_no;
                                    if (k == i)
                                    {
                                        bin_to_add = bin_to_add + j;
                                    }
                                    else
                                    {
                                        bin_to_add = bin_to_add + bins[k];
                                    }
                                }
                                // if (checked_bins.find(bin_to_add) == checked_bins.end())
                                {
                                    auto inserted = bins_to_check_update.insert(bin_to_add);
                                    if (inserted.second)
                                    {
                                        bins_to_check_pair_temp.push_back(std::make_pair(m, bin_to_add));
                                    }
                                    // std::bitset<64> bta_bin(bin_to_add);
                                    // std::cout << bta_bin << "/n";
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (auto itr : (*bins_to_check)[m])
                    {
                        // std::bitset<8*sizeof(itr)> itr_bin(itr);
                        // std::cout << "itr_bin: " << itr_bin << std::endl;
                        // std::cout << "btcu: " << std::endl;
                        // std::cout << "itr: " << itr << "\n";
                        int bins[hyperplane_no];
                        uint64_t temp = itr;
                        for (int i = 0; i < hyperplane_no; i++)
                        {
                            bins[i] = temp % 2;
                            temp = (temp) / 2;
                            // std::cout << "temp: " << temp << "\n";
                        }

                        for (int i = 0; i < hyperplane_no; i++)
                        {

                            // celcing bins to check
                            uint64_t bin_to_add = 0;

                            for (int k = hyperplane_no - 1; k >= 0; k--)
                            {
                                bin_to_add = bin_to_add * 2;
                                if (k == i)
                                {
                                    bin_to_add = bin_to_add + (1 - bins[k]);
                                }
                                else
                                {
                                    bin_to_add = bin_to_add + bins[k];
                                }
                            }
                            // if (checked_bins.find(bin_to_add) == checked_bins.end())
                            {
                                auto inserted = bins_to_check_update.insert(bin_to_add);
                                if (inserted.second)
                                {
                                    bins_to_check_pair_temp.push_back(std::make_pair(m, bin_to_add));
                                }
                                // std::bitset<8*sizeof(bin_to_add)> bta_bin(bin_to_add);
                                // std::cout << bta_bin << std::endl;
                            }
                            // std::cout <<"bin_to_add: " << bin_to_add << "\n";
                        }
                        // std::cout << " itr loop\n";
                    }
                }
                (*bins_to_check)[m] = bins_to_check_update;
            }
            bins_to_check_pair->swap(bins_to_check_pair_temp);
        }

        void render_weighted_hamming(std::vector<std::pair<int, uint64_t>> *bins_to_check_pair, std::vector<std::set<uint64_t>> *bins_to_check, std::vector<std::set<uint64_t>> *checked_bins)
        {

            std::vector<std::pair<int, uint64_t>> bins_to_check_pair_temp;
            std::vector<std::pair<dist_t, std::pair<int, uint64_t>>> bins_to_check_pair_temp_order;
            bins_to_check_pair_temp.reserve(total_table_number * pq_table_number * centroid_no);
            bins_to_check_pair_temp_order.reserve(total_table_number * pq_table_number * centroid_no);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair_temp_hyperplane;
            bins_to_check_pair_temp_hyperplane.reserve(hyperplane_table_number * pq_table_number * centroid_no);
            for (int m = 0; m < total_table_number; m++)
            {
                std::set<uint64_t> bins_to_check_update;
                if (m < pq_table_number)
                {
                    for (auto itr : (*bins_to_check)[m])
                    {
                        // std::bitset<64> itr_bin(itr);
                        // std::cout << "itr_bin: " << itr_bin << "\n";
                        // std::cout << "btcu: \n";
                        int bins[subdivision];
                        uint64_t temp = itr;
                        for (int i = 0; i < subdivision; i++)
                        {
                            bins[i] = temp % centroid_no;
                            temp = (temp) / centroid_no;
                        }

                        for (int i = 0; i < subdivision; i++)
                        {
                            for (int j = 0; j < centroid_no; j++)
                            {
                                // celcing bins to check
                                uint64_t bin_to_add = 0;

                                for (int k = subdivision - 1; k >= 0; k--)
                                {
                                    bin_to_add = bin_to_add * centroid_no;
                                    if (k == i)
                                    {
                                        bin_to_add = bin_to_add + j;
                                    }
                                    else
                                    {
                                        bin_to_add = bin_to_add + bins[k];
                                    }
                                }
                                // if (checked_bins.find(bin_to_add) == checked_bins.end())
                                {
                                    auto inserted = bins_to_check_update.insert(bin_to_add);
                                    if (inserted.second)
                                    {
                                        // std::cout << "bins[i]: " << bins[i] << " j: " << j << " pq_dist_storage[m][i][bins[i]][j]: " << pq_dist_storage[m][i][bins[i]][j] <<"\n";
                                        bins_to_check_pair_temp_order.push_back(std::make_pair(pq_dist_storage[m][i][bins[i]][j], std::make_pair(m, bin_to_add)));
                                    }
                                    // std::bitset<64> bta_bin(bin_to_add);
                                    // std::cout << bta_bin << "/n";
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (auto itr : (*bins_to_check)[m])
                    {
                        // std::bitset<8*sizeof(itr)> itr_bin(itr);
                        // std::cout << "itr_bin: " << itr_bin << std::endl;
                        // std::cout << "btcu: " << std::endl;
                        // std::cout << "itr: " << itr << "\n";
                        int bins[hyperplane_no];
                        uint64_t temp = itr;
                        for (int i = 0; i < hyperplane_no; i++)
                        {
                            bins[i] = temp % 2;
                            temp = (temp) / 2;
                            // std::cout << "temp: " << temp << "\n";
                        }

                        for (int i = 0; i < hyperplane_no; i++)
                        {

                            // celcing bins to check
                            uint64_t bin_to_add = 0;

                            for (int k = hyperplane_no - 1; k >= 0; k--)
                            {
                                bin_to_add = bin_to_add * 2;
                                if (k == i)
                                {
                                    bin_to_add = bin_to_add + (1 - bins[k]);
                                }
                                else
                                {
                                    bin_to_add = bin_to_add + bins[k];
                                }
                            }
                            // if (checked_bins.find(bin_to_add) == checked_bins.end())
                            {
                                auto inserted = bins_to_check_update.insert(bin_to_add);
                                if (inserted.second)
                                {
                                    bins_to_check_pair_temp_hyperplane.push_back(std::make_pair(m, bin_to_add));
                                }
                                // std::bitset<8*sizeof(bin_to_add)> bta_bin(bin_to_add);
                                // std::cout << bta_bin << std::endl;
                            }
                            // std::cout <<"bin_to_add: " << bin_to_add << "\n";
                        }
                        // std::cout << " itr loop\n";
                    }
                }
                (*bins_to_check)[m] = bins_to_check_update;
            }

            std::sort(bins_to_check_pair_temp_order.begin(), bins_to_check_pair_temp_order.end());

            for (auto itr : bins_to_check_pair_temp_order)
            {
                bins_to_check_pair_temp.push_back(itr.second);
                // std::cout << itr.first << " " << itr.second.first << " " << itr.second.second << "\n";
            }
            bins_to_check_pair_temp.insert(bins_to_check_pair_temp.end(), bins_to_check_pair_temp_hyperplane.begin(), bins_to_check_pair_temp_hyperplane.end());

            bins_to_check_pair->swap(bins_to_check_pair_temp);
        }

        void initializeGraph_Point(tableint query, size_t k, size_t min_count = 1, bool pq_dist_applied = false, bool cut_off_applied = false, float cut_off = 1.0)
        { 
            
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>> result;
           
            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair;
          
            std::unordered_set<tableint> checked_ids;

            std::vector<uint64_t> q_bin = getBin(query);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
                bins_to_check_pair.push_back(std::make_pair(i, q_bin[i]));
            
            }

            std::vector<std::set<uint64_t>> checked_bins(total_table_number);
        
            while (top_candidates.size() < min_count)
            {
          
                for (auto itr : bins_to_check_pair)
                {


                    tableint internalId = bin_lookup_[itr.first][itr.second];

                    while (internalId < NO_COND)
                    {
                     
                      //  point_tracker++;
                        bool dup_check = checked_ids.insert(internalId).second;
                       
                        if (dup_check)
                        {
                            dist_t calc_dist = fstdistfunc_(getDataByInternalId(query), getDataByInternalId(internalId), dist_func_param_);
                            top_candidates.emplace(-calc_dist, internalId); 
                          
                        }

                        internalId = *get_link_next(internalId, itr.first);
                       
                    }

                    
                        
                }

          
                if (top_candidates.size() >= min_count)
                {


                    break;
                }

                for (int m = 0; m < total_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
               
                {

                    if (pq_dist_applied)
                    {
             
                        render_weighted_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                    else
                    {

  
                        render_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                }

            }
       

            while (result.size() < k && !(top_candidates.size() == 0))
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
          
                if (rez.second != query)
                {
                    result.push(std::pair<dist_t, tableint>(-rez.first, rez.second));
                }
                top_candidates.pop();
            }
            

            mutuallyConnectNewElement(query, result);
        }

        std::vector<tableint> getConnections(tableint internalId)
        {
            unsigned int *data = get_graph_list(internalId);
            int size = getListCount(data);
            std::vector<tableint> result(size);
            tableint *ll = (tableint *)(data + 1);
            memcpy(result.data(), ll, size * sizeof(tableint));
            return result;
        }

        void initializeGraph_PointRandom(tableint query, size_t k, int64_t cutoff_rand, Xorshift128Plus rng, size_t min_count = 1, bool pq_dist_applied = false, bool cut_off_applied = false, float cut_off = 1.0)
        { 
  

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>> result;
         
            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair;
          
            std::unordered_set<tableint> checked_ids;

            std::vector<uint64_t> q_bin = getBin(query);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
                bins_to_check_pair.push_back(std::make_pair(i, q_bin[i]));
              
            }

        
            std::vector<std::set<uint64_t>> checked_bins(total_table_number);
        
            while (top_candidates.size() < min_count)
            {
               
                for (auto itr : bins_to_check_pair)
                {
                 

                    tableint internalId = bin_lookup_[itr.first][itr.second];
               
                    while (internalId < NO_COND)
                    {
                        if (cutoff_rand < rng.next())
                        {
                            continue;
                        }
                     
                  
                        bool dup_check = checked_ids.insert(internalId).second;
                       
                        if (dup_check)
                        {
                            dist_t calc_dist = fstdistfunc_(getDataByInternalId(query), getDataByInternalId(internalId), dist_func_param_);
                            top_candidates.emplace(-calc_dist, internalId); 
        
                        }

                        internalId = *get_link_next(internalId, itr.first);
                      
                    }

  
                        
                }

                if (top_candidates.size() >= min_count)
                {

                   
                    break;
                }

                for (int m = 0; m < total_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
                {

                    if (pq_dist_applied)
                    {
                        render_weighted_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                    else
                    {

                        render_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                }

            }

            while (result.size() < k && !(top_candidates.size() == 0))
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                if (rez.second != query)
                {
                    result.push(std::pair<dist_t, tableint>(-rez.first, rez.second));
                }
                top_candidates.pop();
            }

            mutuallyConnectNewElement(query, result);
 
        }

        void initializeGraph(bool random_applied = false, float chance = 0.1)
        {
            int count = 0;
            if (random_applied)
            {
                int64_t seed = std::time(0);
                Xorshift128Plus rng(seed);
                uint64_t cutoff_rand = (uint64_t)(chance * UINT64_MAX);
#pragma omp parallel for
                for (size_t i = 0; i < cur_element_count; i++)
                {
      
                    initializeGraph_PointRandom(i, maxM0_, cutoff_rand, rng);
                    
                    {
                        count++;
                    }
                    if (count % 10000 == 0)
                    {
                        std::cout << count << "\n";
                    }
                    
                }
            }
            else
            {

#pragma omp parallel for
                for (size_t i = 0; i < cur_element_count; i++)
                {
                    initializeGraph_Point(i, maxM0_);
                    // #pragma omp critical
                    {
                        count++;
                    }
                    if (count % 10000 == 0)
                    {
                        std::cout << count << "\n";
                    }
                
                }
            }
        }

        template <bool has_deletions>
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayerSTGraph(tableint ep_id, const void *data_point, size_t ef)
        { 
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

            dist_t lowerBound;
            if (!has_deletions || !isMarkedDeleted(ep_id))
            {
                dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
                lowerBound = dist;
                top_candidates.emplace(dist, ep_id);
                candidate_set.emplace(-dist, ep_id);
            }
            else
            {
                lowerBound = std::numeric_limits<dist_t>::max();
                candidate_set.emplace(-lowerBound, ep_id);
            }

            visited_array[ep_id] = visited_array_tag;

            while (!candidate_set.empty())
            {

                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();

                if ((-current_node_pair.first) > lowerBound && (top_candidates.size() == ef || has_deletions == false))
                {
                    break;
                }
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;

                int *data = (int *)get_graph_list(current_node_id);
                size_t size = getListCount((linklistsizeint *)data);
                //              

                for (size_t j = 1; j <= size; j++)
                {
                    int candidate_id = *(data + j);
                    if (!(visited_array[candidate_id] == visited_array_tag))
                    {

                        visited_array[candidate_id] = visited_array_tag;

                        char *currObj1 = (getDataByInternalId(candidate_id));

                        dist_t dist = fstdistfunc_(data_point, currObj1, dist_func_param_);

                        if (top_candidates.size() < ef || lowerBound > dist)
                        {
                            candidate_set.emplace(-dist, candidate_id);

                            if (!has_deletions || !isMarkedDeleted(candidate_id))
                                top_candidates.emplace(dist, candidate_id);

                            if (top_candidates.size() > ef)
                                top_candidates.pop();

                            if (!top_candidates.empty())
                                lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }

            visited_list_pool_->releaseVisitedList(vl);
            return top_candidates;
        }

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn_Graph(const void *query_data, size_t k, size_t ef, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, bool pq_dist_applied, bool cut_off_applied = false, float cut_off = 1.0)
        { 
            if (ef == 0)
            {
                ef = maxM0_;
            }
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            if (cur_element_count == 0)
                return result;
            tableint enterpoint_node_ = 0;
           
            tableint currObj = enterpoint_node_;

        

            dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(enterpoint_node_), dist_func_param_);

            bool changed = true;
            while (changed)
            {
                changed = false;

                unsigned int *data;

                data = (unsigned int *)get_graph_list(currObj);

                int size = getListCount(data);
               

                tableint *datal = (tableint *)(data + 1);
                for (int i = 0; i < size; i++)
                {
                    tableint cand = datal[i];
                    if (cand < 0 || cand > max_elements_)
                    {
                        throw std::runtime_error("cand error");
                    }
                    dist_t d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_);

                    if (d < curdist)
                    {
                        curdist = d;
                        currObj = cand;
                        changed = true;
                    }
                }
            }

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            if (num_deleted_)
            {
                top_candidates = searchBaseLayerSTGraph<true>(
                    currObj, query_data, std::max(ef, k));
            }
            else
            {
                top_candidates = searchBaseLayerSTGraph<false>(
                    currObj, query_data, std::max(ef, k));
            }

            // std::cout << "wowee\n";
            while (top_candidates.size() > k)
            {

                // std::pair<dist_t, tableint> rez = top_candidates.top();
                // std::cout << "tableinr: " << rez.second << "dist_t: " << rez.first << "\n";
                top_candidates.pop();
            }

            while (top_candidates.size() > 0)
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
                top_candidates.pop();
            }
            return result;
        };

        tableint get_enterpoint(const void *query_data)
        { 
            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair;
            
            std::unordered_set<tableint> checked_ids;

            std::vector<uint64_t> q_bin = calculateBin(query_data);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
                bins_to_check_pair.push_back(std::make_pair(i, q_bin[i]));
               
            }

            
            std::vector<std::set<uint64_t>> checked_bins(total_table_number);
            while (true)
            {
           

                for (auto itr : bins_to_check_pair)
                {
                    int point_tracker = 0;
                    int count = 0;
                    if(bin_lookup_[itr.first].find(itr.second) != bin_lookup_[itr.first].end()){
                        return bin_lookup_[itr.first][itr.second];}
                }

                
                for (int m = 0; m < total_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
                {

                    render_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                }
            }
        }

        std::priority_queue<std::pair<dist_t, labeltype>> searchKnn_Hamming(const void *query_data, size_t k, size_t min_count, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, bool pq_dist_applied, bool cut_off_applied = false, float cut_off = 1.0)
        { 
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            if (cur_element_count == 0)
                return result;

     
            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair;
           
            std::unordered_set<tableint> checked_ids;

            std::vector<uint64_t> q_bin = calculateBin(query_data);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
                bins_to_check_pair.push_back(std::make_pair(i, q_bin[i]));
              
            }
            std::vector<std::set<uint64_t>> checked_bins(total_table_number);

            while (top_candidates.size() < min_count)
            {
                for (auto itr : bins_to_check_pair)
                {
                    int point_tracker = 0;
                    checked_bins_no++;
                    int count = 0;

                           
                    tableint internalId = END_OF_BIN;
                    if (bin_lookup_[itr.first].find(itr.second) != bin_lookup_[itr.first].end())
                    {
                        internalId = bin_lookup_[itr.first][itr.second];
                    }

                    while (internalId != END_OF_BIN)
                    {
                        {                           
                           
                            bool dup_check;
                            if (true)
                            {
                                dup_check = checked_ids.insert(internalId).second;
                            }
                            else
                            {
                                dup_check = true;
                            }
                            if (dup_check)
                            {
                                
                                dist_t calc_dist = fstdistfunc_(query_data, getDataByInternalId(internalId), dist_func_param_);
                                top_candidates.emplace(calc_dist, internalId);
                            }


                     }
                        internalId = *get_link_next(internalId, itr.first);
                    }




                }


                if (top_candidates.size() >= min_count)
                {
                    break;
                }

                for (int m = 0; m < total_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
                {

                    if (pq_dist_applied)
                    {
                        render_weighted_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                    else
                    {

                        render_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                }
                for (auto itr : bins_to_check_pair)
                {
                }
            }
            while (top_candidates.size() > k)
            {

                top_candidates.pop();
            }

            while (result.size() < k && !(top_candidates.size() == 0))
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
                top_candidates.pop();
            }
            return result;
        }

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn_ST(const void *query_data, size_t k, size_t ef_, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, bool pq_dist_applied, int cutoff_type = 0, float cutoff = 2.0)
        { 
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            if (cur_element_count == 0)
            {
                return result;
            }

            tableint currObj;

            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair;
            std::unordered_set<tableint> checked_ids;
            std::vector<uint64_t> q_bin = calculateBin(query_data);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
                bins_to_check_pair.push_back(std::make_pair(i, q_bin[i]));

            }
            bool found_viable = false;
            std::vector<std::set<uint64_t>> checked_bins(total_table_number);
            while (!found_viable)
            {

                for (auto itr : bins_to_check_pair)
                {
                    int count = 0;
                    if (bin_lookup_[itr.first].find(itr.second) != bin_lookup_[itr.first].end())
                    {
                        currObj = bin_lookup_[itr.first][itr.second];
                        found_viable = true;
                        break;
                    }
                    if (found_viable)
                    {
                        break;
                    }
                }
                if (found_viable)
                {
                    break;
                }

                for (int m = 0; m < pq_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
                {
                    if (pq_dist_applied)
                    {
                        render_weighted_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                    else
                    {

                        render_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                }
            }


            dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(currObj), dist_func_param_);

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
                  top_candidates = searchBaseLayerST( // done
                currObj, query_data, std::max(ef_, k), checked_points_no, visited_points_no, checked_bins_no, visited_bins_pop, cutoff_type, cutoff);
               while (top_candidates.size() > k)
            {

                      top_candidates.pop();
            }

            while (top_candidates.size() > 0)
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
                top_candidates.pop();
            }
            return result;
        };

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn_ST_inf(const void *query_data, size_t k, size_t ef_, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop, bool pq_dist_applied, int cutoff_type = 1, float cutoff = 1000)
        {
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            if (cur_element_count == 0)
            {
                return result;
            }

            tableint currObj;

            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::vector<std::pair<int, uint64_t>> bins_to_check_pair;
            std::unordered_set<tableint> checked_ids;
            std::vector<uint64_t> q_bin = calculateBin(query_data);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
                bins_to_check_pair.push_back(std::make_pair(i, q_bin[i]));

                         }
            bool found_viable = false;
            std::vector<std::set<uint64_t>> checked_bins(total_table_number);
             while (!found_viable)
            {

                for (auto itr : bins_to_check_pair)
                {
                    int count = 0;
                    auto it = bin_lookup_[itr.first].find(itr.second);
                    if (it != bin_lookup_[itr.first].end())
                    {
                        if (bin_lookup_[itr.first].find(itr.second) != bin_lookup_[itr.first].end())
                        {
                            currObj = bin_lookup_[itr.first][itr.second];
                            found_viable = true;
                            break;
                        }
                    }
                    if (found_viable)
                    {
                        break;
                    }
                }
                if (found_viable)
                {
                    break;
                }

                for (int m = 0; m < pq_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
                {
                    if (pq_dist_applied)
                    {
                        render_weighted_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                    else
                    {

                        render_hamming(&bins_to_check_pair, &bins_to_check, &checked_bins);
                    }
                }
            }


            dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(currObj), dist_func_param_);

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
         
            top_candidates = searchBaseLayerST_inf( 
                currObj, query_data, std::max(ef_, k), checked_points_no, visited_points_no, checked_bins_no, visited_bins_pop, cutoff_type, cutoff);

            while (top_candidates.size() > k)
            {

             
                top_candidates.pop();
            }

            while (top_candidates.size() > 0)
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
                top_candidates.pop();
            }
            return result;
        };

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn_ST_nh(const void *query_data, size_t k, size_t ef_, int &checked_points_no, int &visited_points_no, int &checked_bins_no, std::vector<int> *visited_bins_pop)
        { 
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            if (cur_element_count == 0)
            {
                return result;
            }

            tableint currObj;

            std::vector<std::set<uint64_t>> bins_to_check(total_table_number);
            std::unordered_set<tableint> checked_ids;
            std::vector<uint64_t> q_bin = calculateBin(query_data);
            for (int i = 0; i < total_table_number; i++)
            {
                bins_to_check[i].insert(q_bin[i]);
             
            }
            bool found_viable = false;
            std::vector<std::set<uint64_t>> checked_bins(total_table_number);
       
            while (!found_viable)
            {
                for (int i = 0; i < total_table_number; i++)
                {
                    for (auto itr : bins_to_check[i])
                    {

                          
                            if (bin_lookup_[i].find(itr) != bin_lookup_[i].end())
                            {
                                currObj = bin_lookup_[i][itr];
                                found_viable = true;
                                break;
                            }
                    }
                    if (found_viable)
                    {
                        break;
                    }
                }
                if (found_viable)
                {
                    break;
                }
                else
                {
                    return result;
                }

                for (int m = 0; m < pq_table_number; m++)
                {
                    for (auto itr : (bins_to_check)[m])
                    {
                        (checked_bins)[m].insert(itr);
                    }
                }
            }

           

            dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(currObj), dist_func_param_);

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;

            top_candidates = searchBaseLayerST(
                currObj, query_data, std::max(ef_, k), checked_points_no, visited_points_no, checked_bins_no, visited_bins_pop);

       
            while (top_candidates.size() > k)
            {

           
                top_candidates.pop();
            }

            while (top_candidates.size() > 0)
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
                top_candidates.pop();
            }
            return result;
        };

        void setListCount(linklistsizeint *ptr, unsigned short int size) const
        {
            *((unsigned short int *)(ptr)) = *((unsigned short int *)&size);
        }

        void setBin(tableint internalId, std::vector<uint64_t> bin_val)
        {
            for (int i = 0; i < total_table_number; i++)
            {
                unsigned char *bin_cur = (get_binlist(internalId, i));
                memcpy(bin_cur, &(bin_val[i]), sizeof(u_int64_t)); // done
                add_bin_lookup(internalId, bin_val[i], i);
            }
        }

        void setBin(tableint internalId, uint64_t bin_val, int pq_table_no)
        {
            unsigned char *bin_cur = (get_binlist(internalId, pq_table_no)); // done
            memcpy(bin_cur, &bin_val, sizeof(u_int64_t));                    // done

            add_bin_lookup(internalId, bin_val, pq_table_no);
        }

        const int traverseThroughTables(){
            size_t totalCount = 0;
            for (size_t table_no = 0; table_no < bin_lookup_.size(); ++table_no) {
                const auto& umap = bin_lookup_[table_no];
                for (const auto& kvp : umap) {
                    tableint internal_id = kvp.second; // Assuming the value in the map is the starting element of the linked list
                    while (internal_id < UINT_MAX - 3) {
                        totalCount++;
                        internal_id = *get_link_next(internal_id, table_no);
                    }
                }
            }
            return totalCount;
        }

        void getNeighborsByHeuristic2(
            std::priority_queue<std::pair<dist_t, tableint>> &top_candidates,
            const size_t M)
        {
            if (top_candidates.size() < M)
            {
                return;
            }

            std::priority_queue<std::pair<dist_t, tableint>> queue_closest;
            std::vector<std::pair<dist_t, tableint>> return_list;
            while (top_candidates.size() > 0)
            {
                queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
                top_candidates.pop();
            }

            while (queue_closest.size())
            {
                if (return_list.size() >= M)
                    break;
                std::pair<dist_t, tableint> curent_pair = queue_closest.top();
                dist_t dist_to_query = -curent_pair.first;
                queue_closest.pop();
                bool good = true;

                for (std::pair<dist_t, tableint> second_pair : return_list)
                {
                    dist_t curdist =
                        fstdistfunc_(getDataByInternalId(second_pair.second),
                                     getDataByInternalId(curent_pair.second),
                                     dist_func_param_);
                    ;
                    if (curdist < dist_to_query)
                    {
                        good = false;
                        break;
                    }
                }
                if (good)
                {
                    return_list.push_back(curent_pair);
                }
            }

            for (std::pair<dist_t, tableint> curent_pair : return_list)
            {
                top_candidates.emplace(-curent_pair.first, curent_pair.second);
            }
        }

        void mutuallyConnectNewElement(tableint cur_c,
                                       std::priority_queue<std::pair<dist_t, tableint>> &top_candidates)
        {

            size_t Mcurmax = maxM0_;

            std::vector<tableint> selectedNeighbors;
            selectedNeighbors.reserve(Mcurmax);
            while (top_candidates.size() > 0 && Mcurmax > selectedNeighbors.size())
            {
                selectedNeighbors.push_back(top_candidates.top().second);
                top_candidates.pop();
            }

            {
                linklistsizeint *ll_cur;
                ll_cur = get_graph_list(cur_c);
                setListCount(ll_cur, selectedNeighbors.size());
                tableint *data = (tableint *)(ll_cur + 1);
                for (size_t idx = 0; idx < selectedNeighbors.size(); idx++)
                {

                    data[idx] = selectedNeighbors[idx];
                }
            }
                   }

        template <bool setbin = true>
        std::vector<uint64_t> updatePoint(const void *dataPoint, tableint internalId)
        { //
            throw std::runtime_error("updateing not currently vaild");
            std::vector<uint64_t> newbin;

            memcpy(getDataByInternalId(internalId), dataPoint, data_size_);
            unsigned char *bin_cur = (get_binlist(internalId)) + 2;

            for (int i = 0; i < total_table_number; i++)
            {

                delete_bin_lookup(internalId, i);
            }

            if (setbin)
            {
                newbin = calculateBin(dataPoint);
                setBin(internalId, newbin);
                // std::cout << "setbin\n";

                if (!isMarkedBinSet(internalId))
                {
                    markBinSetInternal(internalId);
                };
            }
            else
            {

                if (isMarkedBinSet(internalId))
                {
                    unmarkBinSetInternal(internalId);
                };
            }

            return newbin;
        };

        void checkBinsLength(int pq_table_no = 0)
        {
            int total = 0;
            int bin_max = 0;
            for (int pq_table_no = 0; pq_table_no < pq_table_number; pq_table_no++)
            {
                if (myPow(centroid_no, subdivision) < myPow(2, hyperplane_no))
                {
                    bin_max = myPow(2, hyperplane_no);
                }
                else
                {
                    bin_max = myPow(centroid_no, subdivision);
                }
                for (int i = 0; i < bin_max; i++)
                {
                    int count = 0;


                    tableint internalId = END_OF_BIN;
                    if (bin_lookup_[pq_table_no].find(i) != bin_lookup_[pq_table_no].end())
                    {   
                        internalId = bin_lookup_[pq_table_no][i];
                    }

                    while (internalId != END_OF_BIN)
                    {   
                        {
                        count++;




                     }
                        internalId = *get_link_next(internalId, pq_table_no);
                    }


                    total = total + count;
                }
            }
            std::cout << "bin total:: " << total << "\n";
        }

        void initializeBins()
        {
            int count = 0;
            cutoff_saved = 0.0;

#pragma omp parallel for
            for (int internalId = 0; internalId < cur_element_count; internalId++)
            {
                if (!isMarkedBinSet(internalId))
                {
                    std::vector<uint64_t> newbin = calculateBin(getDataByInternalId(internalId));

                    unsigned char *bin_cur = (get_binlist(internalId)) + 2; // done
                    for (int i = 0; i < total_table_number; i++)
                    {
                        memcpy(bin_cur + i * sizeof(u_int64_t), &(newbin[i]), sizeof(u_int64_t));
                        *get_link_next(internalId, i) = DISCONNECTED;
                        *get_link_prev(internalId, i) = DISCONNECTED;

                        add_bin_lookup(internalId, newbin[i], i);
                    }

                    markBinSetInternal(internalId);
                }
            }
                   bins_initialized = 1;
        }

        void initializeCentroidList(bool rand = false, std::string metric = "", std::string learn_data = "") // TODO
        {
            std::uniform_real_distribution<float> dist(-1.0, 1.0);
            std::default_random_engine random_generator_;
            std::vector<int> vec(cur_element_count);
            std::iota(std::begin(vec), std::end(vec), 0);
            std::shuffle(std::begin(vec), std::end(vec), random_generator_);
            std::cout << "pq_t_n: " << pq_table_number << ", subdivision: " << subdivision << ", centroid_no: " << centroid_no << "\n";
            if (rand)
            {

                for (int k = 0; k < pq_table_number; k++)
                {
                 
#pragma omp parallel for
                    for (int i = 0; i < subdivision; i++)
                    {
                        for (int j = 0; j < centroid_no; j++)
                        {
                            memcpy(get_centroidlist(k, j, i), getDataByInternalId(vec[k * subdivision * centroid_no + i * centroid_no + j]) + data_size_sub_ * i, data_size_sub_);

                         
                        }
                    }
                }
            }
            // else
            // {

            //     double t0 = elapsed();
            //     if (metric == "euclidean")
            //     {

            //         size_t vecdim_t = vecdim;
            //         size_t nta;
            //         float *xt = fvecs_read(learn_data.c_str(), &vecdim_t, &nta);
            //         std::cout << "nta: " << nta << "\n";
            //         // #pragma omp parallel
            //         for (int i = 0; i < pq_table_number; i++)
            //         {
            //             size_t nt = nta / pq_table_number;
            //             // nt = 2;
            //             std::cout << "i: " << i << "\n";
            //             std::cout << "nt: " << nt << "\n";
            //             std::cout << "vecdim: " << vecdim << "\n";

            //             int ksub = std::ceil(log2(centroid_no));
            //             std::cout << "ksub: " << ksub << "\n";
            //             faiss::IndexPQ index(vecdim, subdivision, ksub, faiss::METRIC_L2);

            //             printf("[%.3f s] Loading train set\n", elapsed() - t0);

            //             float *xt_c = xt + (nt * i) * (vecdim_t);
            //             /*std::cout << "xt:" << xt[0]  << "\n";
            //             std::cout << " addr: " << xt << "\n";
            //             std::cout << "xt_c:" << xt_c[0] << "\n";
            //             std::cout << "xt[nt]:" << xt[nt*i*vecdim_t] << "\n";
            //             std::cout << " addr: " << xt_c << "\n";*/

            //             printf("[%.3f s] Preparing index \"%ld\" d=%ld\n", elapsed() - t0, i, vecdim);

            //             std::cout << "centroids_size: " << index.pq.centroids.size() << "\n";
            //             /*for (int l = 0; l < index.pq.centroids.size(); l++)
            //             {

            //                 std::cout << index.pq.centroids[l] << "\t";
            //             }*/
            //             std::cout << "\n";
            //             printf("[%.3f s] Training on %ld vectors\n", elapsed() - t0, nt);
            //             index.verbose = true;
            //             // std::cout << "asasfag\n";
            //             index.train(nt, xt_c);
            //             // std::cout << "aaa\n";
            //             //  float* array_centroids = &index.pq.centroids[0];
            //             // std::cout << "check prev\n";
            //             for (int j = 0; j < subdivision; j++)
            //             {
            //                 for (int k = 0; k < centroid_no; k++)
            //                 {
            //                     memcpy(get_centroidlist(i, k, j), index.pq.get_centroids(j, k), data_size_sub_);
            //                     /*for (int m = 0; m < data_size_sub_ / sizeof(float); m++)
            //                     {
            //                         char temp[4];
            //                         for (int l = 0; l < 4; l++)
            //                         {
            //                             temp[l] = ((char *)get_centroidlist(j, k, i)+sizeof(float)*m)[l];
            //                         }
            //                         float *temp_float = (float *)temp;
            //                         std::cout << i << " " << j << " " << k << " " << *temp_float << " " << (index.pq.get_centroids(j, k))[m] << "\n";
            //                     }*/
            //                 }
            //             }

            //             // float *lol = index.pq.get_centroids(1, 1);
            //             // float wow = (float)*lol;
            //             // std::cout << "wow" << wow << "\n";
            //             // std::cout << "centroids_size: " << index.pq.centroids.size() << "\n";
            //             // std::cout << "M: " << index.pq.M << " ksub: " << index.pq.ksub << " dsub: " << index.pq.dsub << "\n";
            //             /*for (int l = 0; l < index.pq.centroids.size(); l++)
            //             {

            //                 std::cout << l + 1 << ": " << index.pq.centroids[l] << "\t";
            //             }
            //             std::cout << "\n";
            //             for (int l = 0; l < index.pq.transposed_centroids.size(); l++)
            //             {

            //                 std::cout << l + 1 << ": " << index.pq.transposed_centroids[l] << "\t";
            //             }*/
            //         }

            //         delete[] xt;
            //     }
            //     else if (metric == "cosine" || metric == "angular")
            //     {

            //         size_t vecdim_t = vecdim;
            //         size_t nta;
            //         float *xt = fvecs_read(learn_data.c_str(), &vecdim_t, &nta);
            //         std::cout << "nta: " << nta << "\n";
            //         // #pragma omp parallel
            //         for (int i = 0; i < pq_table_number; i++)
            //         {
            //             size_t nt = nta / pq_table_number;
            //             // nt = 2;
            //             std::cout << "i: " << i << "\n";
            //             std::cout << "nt: " << nt << "\n";

            //             int ksub = std::ceil(log2(centroid_no));
            //             std::cout << "ksub: " << ksub << "\n";
            //             faiss::IndexPQ index(vecdim, subdivision, ksub, faiss::METRIC_INNER_PRODUCT);

            //             printf("[%.3f s] Loading train set\n", elapsed() - t0);

            //             float *xt_c = xt + (nt * i) * (vecdim_t);
            //             /*std::cout << "xt:" << xt[0]  << "\n";
            //             std::cout << " addr: " << xt << "\n";
            //             std::cout << "xt_c:" << xt_c[0] << "\n";
            //             std::cout << "xt[nt]:" << xt[nt*i*vecdim_t] << "\n";
            //             std::cout << " addr: " << xt_c << "\n";*/

            //             printf("[%.3f s] Preparing index \"%ld\" d=%ld\n", elapsed() - t0, i, vecdim);

            //             std::cout << "centroids_size: " << index.pq.centroids.size() << "\n";
            //             /*for (int l = 0; l < index.pq.centroids.size(); l++)
            //             {

            //                 std::cout << index.pq.centroids[l] << "\t";
            //             }*/
            //             std::cout << "\n";
            //             printf("[%.3f s] Training on %ld vectors\n", elapsed() - t0, nt);
            //             index.verbose = true;
            //             index.train(nt, xt_c);
            //             float *array_centroids = &index.pq.centroids[0];
            //             std::cout << "check prev\n";
            //             for (int j = 0; j < subdivision; j++)
            //             {
            //                 for (int k = 0; k < centroid_no; k++)
            //                 {
            //                     memcpy(get_centroidlist(i, k, j), index.pq.get_centroids(j, k), data_size_sub_);
            //                     /*for (int m = 0; m < data_size_sub_ / sizeof(float); m++)
            //                     {
            //                         char temp[4];
            //                         for (int l = 0; l < 4; l++)
            //                         {
            //                             temp[l] = ((char *)get_centroidlist(j, k, i)+sizeof(float)*m)[l];
            //                         }
            //                         float *temp_float = (float *)temp;
            //                         std::cout << i << " " << j << " " << k << " " << *temp_float << " " << (index.pq.get_centroids(j, k))[m] << "\n";
            //                     }*/
            //                 }
            //             }

            //             float *lol = index.pq.get_centroids(1, 1);
            //             float wow = (float)*lol;
            //             // std::cout << "centroids_size: " << index.pq.centroids.size() << "\n";
            //             // std::cout << "M: " << index.pq.M << " ksub: " << index.pq.ksub << " dsub: " << index.pq.dsub << "\n";
            //             /*for (int l = 0; l < index.pq.centroids.size(); l++)
            //             {

            //                 std::cout << l + 1 << ": " << index.pq.centroids[l] << "\t";
            //             }
            //             std::cout << "\n";
            //             for (int l = 0; l < index.pq.transposed_centroids.size(); l++)
            //             {

            //                 std::cout << l + 1 << ": " << index.pq.transposed_centroids[l] << "\t";
            //             }*/
            //         }

            //         delete[] xt;
            //     }
            // }
            // // int aa;
            // // std::cin >> aa;

            for (int i = pq_table_number; i < total_table_number; i++)
            {
               
                for (int j = 0; j < hyperplane_no; j++)
                {
                    std::vector<float> hyperplane(vecdim);
                    for (int k = 0; k < vecdim; k++)
                    {
                        hyperplane[k] = dist(random_generator_);
                    }
                    memcpy(get_centroidlist(i, j, 0), hyperplane.data(), data_size_);
                }
            }
        }

        void initializeCentroidDistanceStorage()
        {
#pragma omp parallel for collapse(2)
            for (int pq_t_n = 0; pq_t_n < pq_table_number; pq_t_n++)
            {
                for (int subdiv = 0; subdiv < subdivision; subdiv++)
                {
                    for (int a = 0; a < centroid_no; a++)
                    {
                        for (int b = a + 1; b < centroid_no; b++)
                        {
                            pq_dist_storage[pq_t_n][subdiv][a][b] = (dist_t)fstdistfunc_sub_(get_centroidlist(pq_t_n, a, subdiv), get_centroidlist(pq_t_n, b, subdiv), dist_func_param_sub_);
                            pq_dist_storage[pq_t_n][subdiv][b][a] = pq_dist_storage[pq_t_n][subdiv][a][b];
                        }
                    }
                    for (int a = 0; a < centroid_no; a++)
                    {
                        for (int b = 0; b < centroid_no; b++)
                        {
                        }
                    }
                }
            }
        }

        template <bool setbin = false>
        tableint addPoint(const void *data_point, labeltype label, bool cut_off_applied = false, float cut_off = 1.0)
        {

            tableint cur_c = 0;
            {
                std::unique_lock<std::mutex> templock_curr(cur_element_count_guard_);
                auto search = label_lookup_.find(label);
                if (search != label_lookup_.end())
                {
                    tableint existingInternalId = search->second;

                    if (isMarkedDeleted(existingInternalId))
                    {
                        unmarkDeletedInternal(existingInternalId);
                    }
                    // done
                    if (setbin)
                    {
                        updatePoint<true>(data_point, existingInternalId);
                    }
                    else
                    {
                        updatePoint<false>(data_point, existingInternalId);
                    }

                    return existingInternalId;
                }

                if (cur_element_count >= max_elements_)
                {
                    throw std::runtime_error("The number of elements exceeds the specified limit");
                };

                cur_c = cur_element_count;
                cur_element_count++;
                label_lookup_[label] = cur_c;
            }

            insert_time[cur_c] = 0;
            delete_time_indiv[cur_c] = 0;

            //
            memset(data_level0_memory_ + cur_c * size_data_per_element_, 0, size_data_per_element_);

            memcpy(getExternalLabeLp(cur_c), &label, sizeof(labeltype));
            memcpy(getDataByInternalId(cur_c), data_point, data_size_);

            if (bins_initialized == 1)
            {
                if (!isMarkedBinSet(cur_c))
                {
                    std::vector<uint64_t> newbin = calculateBin(getDataByInternalId(cur_c));
                    unsigned char *bin_cur = (get_binlist(cur_c)) + 2; // done
                    for (int i = 0; i < total_table_number; i++)
                    {

                        memcpy(bin_cur + i * sizeof(u_int64_t), &(newbin[i]), sizeof(u_int64_t));
                        *get_link_next(cur_c, i) = DISCONNECTED;
                        *get_link_prev(cur_c, i) = DISCONNECTED;
                    }

                    {
#pragma omp parallel for
                        for (int i = 0; i < total_table_number; i++)
                        {
                            // count = count++;

                            // std::cout << "addingf bin lookup\n";
                            add_bin_lookup(cur_c, newbin[i], i);
                            // std::cout << "finished adding bin lookup\n";
                        }
                    }
                    // std::cout << "marking bin set\n";
                    markBinSetInternal(cur_c);
                    // std::cout << "finshed marking bin set\n";
                }
            }

            return cur_c;
        };

        unsigned int addPoint(const void *datapoint, labeltype label, bool satis)
        {
            tableint temp = addPoint(datapoint, label);
            return temp;
        }

        int get_cur_element_count()
        {
            return cur_element_count;
        };
    };

};
