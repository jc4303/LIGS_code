#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include "ligslib/hnswlib.h"
#include <cmath>
#include <limits>
#include <iomanip>
#include <unordered_set>
#include <sstream>
#include <string>
#include <math.h>
#include <sys/stat.h>
#include <bits/stdc++.h>
#include <filesystem>
#include <random>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>   
#include <pybind11/numpy.h> 
#include <utility>
#include <vector>
#include <queue>
#include <tuple>

#define BOOTSTRAP 100000
#define PARTITION 100001
#define DELETEADD 100002
#define STANDARD 100003
#define BATINSSEQ 100004
#define BATINSBUN 100005
#define BIGSAVE 100006
#define ADDITIONAL 100007
#define INITIALIZE 100008
#define PQSTANDARD 100009
#define PQOVERLAP 100010
#define PQDIST 100011
#define PARTIAL 100012
#define SELECT 100013
#define CLOSEPOINT 100014
#define UPDATE 100015
#define REPLACE 100016
#define REPLACEORIG 100017
#define REPLACENORECON 100018
#define REPLACEZERO 100019
#define STANDARDZERO 100020
#define NOUPDATE 100021
#define ADDPOINT 100022
#define REPLACEADDPOINT 100023
#define ADDSIMPLE 100024
#define ADDORIG 100025
// #define STANDARD 358383
// #define BOOTSTRAPP 142555
// #define PARTITIONP 743366
// #define ADDCOMBI 536378

#define FIFO 1000
#define LIFO 1001
#define RANDOM 1002
#define CORREL 1003
#define CORRELCH 1004

#define NONE 100
#define ALL 101
#define CHECKPOINT 102

#define LAZY 10000
#define MARKDEL 10001
#define MDUPDATE 10002
#define MDUPPRUNE 10003
#define MDUPDATEND 10004
#define MDUPPRUNEND 10005
#define REBUILD 10006

#define LOADGT 200000
#define SAVEGT 200001
#define DEFAULTGT 200002

using namespace std;
using namespace hnswlib;
namespace py = pybind11;





class StopW
{
    std::chrono::steady_clock::time_point time_begin;

public:
    StopW()
    {
        time_begin = std::chrono::steady_clock::now();
    }

    float getElapsedTimeMicro()
    {
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    void reset()
    {
        time_begin = std::chrono::steady_clock::now();
    }
};

//############################# From HNSWlib at https://github.com/nmslib/hnswlib ##################################//
/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))

#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
static size_t getPeakRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L; /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t)0L; /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L; /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
static size_t getCurrentRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return (size_t)0L; /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE *fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        return (size_t)0L; /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1)
    {
        fclose(fp);
        return (size_t)0L; /* Can't read? */
    }
    fclose(fp);
    return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L; /* Unsupported. */
#endif
}

//############################# From HNSWlib at https://github.com/nmslib/hnswlib ##################################//
static void
get_gt(unsigned int *massQA, unsigned char *massQ, size_t vecsize, size_t qsize, SpaceInterface<float> *l2space,
       size_t vecdim, vector<std::priority_queue<std::pair<float, labeltype>>> &answers, size_t k)
{

    (vector<std::priority_queue<std::pair<float, labeltype>>>(qsize)).swap(answers);
    DISTFUNC<float> fstdistfunc_ = l2space->get_dist_func();

    for (int i = 0; i < qsize; i++)
    {

        for (int j = 0; j < k; j++)
        {
            answers[i].emplace(0.0f, massQA[100 * i * 4 + j * 4]);

        }

    }
}

void tokenize(std::string const &str, const char delim,
              std::vector<std::string> &out)
{
    // construct a stream from the string
    std::stringstream ss(str);

    std::string s;
    while (std::getline(ss, s, delim))
    {
        out.push_back(s);
    }
}

static void
calc_gt(unsigned char *massQ, size_t vecsize, size_t qsize, PQbinDoubleGraph<float> *appr_alg, SpaceInterface<float> *l2space,
        size_t vecdim, vector<std::priority_queue<std::pair<float, labeltype>>> &answers)
{

    (vector<std::priority_queue<std::pair<float, labeltype>>>(qsize)).swap(answers);
    DISTFUNC<float> fstdistfunc_ = l2space->get_dist_func();
#pragma omp parallel for
    for (int j = 0; j < qsize; j++)
    {
        // std::cout << "calcuting " << j << "th query\n";
        float distance = std::numeric_limits<float>::max();
        labeltype NN;
        for (int i = 0; i < vecsize; i++)
        {
            if (appr_alg->isLabel((size_t)i))
            {
                float temp_dist = fstdistfunc_(massQ + vecdim * j * 4, appr_alg->getDataByLabelChar((size_t)i), appr_alg->dist_func_param_);
                // std::cout << "label: " << i << "dist: " << temp_dist << "\n";
                if (temp_dist < distance)
                {
                    distance = temp_dist;
                    NN = i;
                }
            }
        }
        answers[j].emplace(distance, NN);
    }
}

static void normalize_angular(unsigned char *mass, int vecdim)
{
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
        norm_vec[i / 4] = (*temp_float);
        total = total + norm_vec[i / 4] * norm_vec[i / 4];
    }
    float total_sqrt = sqrt(total);
    for (int i = 0; i < vecdim; i++)
    {
        char temp[4];
        norm_vec[i] = norm_vec[i] / total_sqrt;
        memcpy(mass, &(norm_vec[i]), sizeof(float));
        for (int j = 0; j < 4; j++)
        {
            mass[i * 4 + j] = temp[j];
        }
    }
}

int myPow(int x, unsigned int p)
{
    if (p == 0)
        return 1;
    if (p == 1)
        return x;

    int tmp = myPow(x, p / 2);
    if (p % 2 == 0)
        return tmp * tmp;
    else
        return x * tmp * tmp;
}

static std::tuple<float, float>
test_approx(unsigned char *massQ, size_t vecsize, size_t qsize, PQbinDoubleGraph<float> &appr_alg, size_t vecdim,
            vector<std::priority_queue<std::pair<float, labeltype>>> &answers, size_t min_count, size_t k, float &not_deleted, float &deleted, float &checked_points_no, float &visited_points_no, float &checked_bins_no, std::vector<int> *visited_bins_pop, bool hamming, bool similarity, bool knn_cutoff_applied, float knn_cutoff, bool weighted, std::string query_type)
{
    std::vector<float> run_times(qsize, 0.0f);
    size_t correct = 0;
    size_t total = 0;
    float dist_total = 0.0;
    float recall = 0.0;

    if (similarity)
    {
        int total_deleted = 0;
        int total_not_deleted = 0;
#pragma omp parallel for
        for (int i = 0; i < qsize; i++)
        {
            int single_deleted = 0;
            int single_not_deleted = 0;
            int checked_points_no_singular = 0;
            int visited_points_no_singular = 0;
            int checked_bins_no_singular = 0;
            std::priority_queue<std::pair<float, labeltype>> result;
            if (hamming)
            {
                result = appr_alg.searchKnn_Hamming(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
            }
            else
            {
                if (query_type == "LIGS_cutoff")
                {
                    result = appr_alg.searchKnn_ST_inf(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }

                else if (query_type == "graph")
                {
                    result = appr_alg.searchKnn_Graph(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
                else
                {
                    result = appr_alg.searchKnn_ST(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
            }
        std::priority_queue<std::pair<float, labeltype>> gt(answers[i]);
        }
    }
    else
    {
        int total_deleted = 0;
        int total_not_deleted = 0;

    
#pragma omp parallel for // done
        for (int i = 0; i < qsize; i++)
        {
            StopW indiv_search_time = StopW();

            int single_deleted = 0;
            int single_not_deleted = 0;
            int checked_points_no_singular = 0;
            int visited_points_no_singular = 0;
            int checked_bins_no_singular = 0;
            std::priority_queue<std::pair<float, labeltype>> result;
            if (hamming)
            {
                result = appr_alg.searchKnn_Hamming(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);

            }
            else
            {

                if (query_type == "LIGS_cutoff")
                {
                    result = appr_alg.searchKnn_ST_inf(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }

                else if (query_type == "graph")
                {
                    result = appr_alg.searchKnn_Graph(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
                else
                {
                    result = appr_alg.searchKnn_ST(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
            }

            std::priority_queue<std::pair<float, labeltype>> gt(answers[i]);

            unordered_set<labeltype> g;
            total += gt.size();

            while (gt.size())
            {

                g.insert(gt.top().second);
                gt.pop();
            }

            while (result.size())
            {
                if (g.find(result.top().second) != g.end())
                {

                    correct++;
                }
                else
                {
                }
                result.pop();
            }
            run_times[i] = indiv_search_time.getElapsedTimeMicro();
        }

        /*std::string filename = "search_times_ligs_" + std::to_string(k)  + "_" + std::to_string(((1.0f * correct / total)*1000.0))  + ".txt";

        // Open the file in write mode
        std::ofstream file(filename);

        // Check if the file opened successfully
        if (!file.is_open()) {
            std::cerr << "Failed to open the file for writing." << std::endl;
            //return 1;
        }

        // Write the data
        for (float time : run_times) {
            file << time << std::endl;
        }

        // Close the file
        file.close();
        std::cout << "Data written to " << filename << std::endl;*/

        recall = 1.0f * correct / total;
    }



    if (similarity)
    {

    StopW full_delete_time = StopW();
    int single_deleted = 0;
    int single_not_deleted = 0;
    int checked_points_no_singular = 0;
    int visited_points_no_singular = 0;
    int checked_bins_no_singular = 0;
    std::priority_queue<std::pair<float, labeltype>> result;
#pragma omp parallel for
        for (int i = 0; i < qsize; i++)
        {
            if (hamming)
            {
                result = appr_alg.searchKnn_Hamming(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
            }
            else
            {
                if (query_type == "LIGS_cutoff")
                {
                    result = appr_alg.searchKnn_ST_inf(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }

                else if (query_type == "graph")
                {
                    result = appr_alg.searchKnn_Graph(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
                else
                {
                    result = appr_alg.searchKnn_ST(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
            }
        }
        float total_search_time = full_delete_time.getElapsedTimeMicro();
        std::cout << "search time: " << total_search_time << "\n";
        return std::make_tuple(1.0f * dist_total / total, total_search_time); //jwc = fix
    }
    else
    {

    StopW full_delete_time = StopW();
    int single_deleted = 0;
    int single_not_deleted = 0;
    int checked_points_no_singular = 0;
    int visited_points_no_singular = 0;
    int checked_bins_no_singular = 0;
    std::priority_queue<std::pair<float, labeltype>> result;
    
#pragma omp parallel for // done
        for (int i = 0; i < qsize; i++)
        {
            std::priority_queue<std::pair<float, labeltype>> result;
            if (hamming)
            {
                result = appr_alg.searchKnn_Hamming(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);

            }
            else
            {

                if (query_type == "LIGS_cutoff")
                {
                    result = appr_alg.searchKnn_ST_inf(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }

                else if (query_type == "graph")
                {
                    result = appr_alg.searchKnn_Graph(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
                else
                {
                    result = appr_alg.searchKnn_ST(massQ + vecdim * i * 4, k, min_count, checked_points_no_singular, visited_points_no_singular, checked_bins_no_singular, visited_bins_pop, weighted);
                }
            }

         
        }
        float total_search_time = full_delete_time.getElapsedTimeMicro();
        std::cout << "search time: " << total_search_time << "\n";
        return std::make_tuple(recall, total_search_time);
    }

}

static std::tuple<std::vector<size_t>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> *
test_vs_recall(unsigned char *massQ, size_t vecsize, size_t qsize, PQbinDoubleGraph<float> &appr_alg, size_t vecdim,
               vector<std::priority_queue<std::pair<float, labeltype>>> &answers, size_t k, bool similarity, bool knn_cutoff_applied, float knn_cutoff, std::string index_file, bool hamming, bool weighted, std::string query_type)
{
    /* for (const auto &pair : appr_alg.bin_lookup_[0]) {
         std::cout << pair.first << ": " << pair.second << std::endl;
     }
         for (const auto &pair : appr_alg.bin_lookup_[1]) {
         std::cout << pair.first << ": " << pair.second << std::endl;
     }&*/
    vector<size_t> min_count_list; // = { 10,10,10,10,10 };
    vector<float> recalls;
    vector<float> qps;
    vector<float> search_deleted;
    vector<float> search_not_deleted;
    vector<float> checked_points_no;
    vector<float> visited_points_no;
    vector<float> checked_bins_no;
    uint64_t max = std::max(myPow(2, appr_alg.hyperplane_no), myPow(appr_alg.centroid_no, appr_alg.subdivision));
    int lim_ = (int)((float)appr_alg.cur_element_count / (double)(max)) / 2.0;
    std::cout << "max: " << max << "\n";
    std::cout << "lim: " << lim_ << "\n";
    std::vector<std::vector<int> *> visited_bins_pop;
    std::vector<int> min_count_to_push = {1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 15, 17, 18, 20, 30, 35, 40, 46, 50, 70, 100, 120, 150, 170, 200, 300, 500};
    if (query_type == "LIGS_cutoff")
    {
        std::cout << "\n"
                  << "using LIGS w/ cutoff \n \n";
    }

    else if (query_type == "graph")
    {
        std::cout << "\n"
                  << "using graph search \n \n";
    }

    else if (query_type == "hamming")
    {
        std::cout << "\n"
                  << "using hamming search \n \n";
    }
    else
    {
        std::cout << "\n"
                  << "using std. LIGS \n \n";
    }

    for (int i = 0; i < min_count_to_push.size(); i++)
    {
        min_count_list.push_back(min_count_to_push[i]);
        search_deleted.push_back(0);
        search_not_deleted.push_back(0);
        checked_points_no.push_back(0);
        visited_points_no.push_back(0);
        checked_bins_no.push_back(0);
        visited_bins_pop.push_back(new std::vector<int>(vecsize));
    }


    std::cout << "hamming: " << hamming << ", weighted: " << weighted << "\n";
    for (int i = 0; i < min_count_list.size(); i++)
    {
        int min_count = min_count_list[i];
        StopW stopw = StopW();
        // std::cout << "check 1 \n";
         auto [recall, total_time ] = test_approx(massQ, vecsize, qsize, appr_alg, vecdim, answers, min_count, k, search_not_deleted[i], search_deleted[i], checked_points_no[i], visited_points_no[i], checked_bins_no[i], visited_bins_pop[i], hamming, similarity, knn_cutoff_applied, knn_cutoff, weighted, query_type);
        float time_us_per_query = total_time / qsize;

        // if(i > 2){
        std::cout << min_count << "\t" << recall << "\t" << time_us_per_query << " us\n";
        // std::cout << ef << "\t" << recall << "\t" << time_us_per_query << " us\n";
        /*if (recall > 1.0) {
            std::cout << recall << "\t" << time_us_per_query << " us\n";
            break;
        }*/
        recalls.push_back(recall);
        qps.push_back(1e+6 / time_us_per_query);
        if ((time_us_per_query) > 20000)
        {
            break;
        }
        //}
    }
    float avg_dist_ep = 0;
    float avg_dist_op = 0;
    float avg_insert_time = 0;
    int no_ = 0;
    for (int i = 0; i < appr_alg.cur_element_count; i++)
    {
        if (!appr_alg.isMarkedDeleted(i))
        {
            avg_insert_time = avg_insert_time + appr_alg.insert_time[i];
            no_++;
        }
    }
    avg_dist_ep = avg_dist_ep / no_;
    avg_dist_op = avg_dist_op / no_;
    avg_insert_time = avg_insert_time / no_;
    std::cout << "avg_dist_ep: " << avg_dist_ep << "\n";
    std::cout << "avg_dist_op: " << avg_dist_op << "\n";
    std::cout << "avg_insert_time: " << avg_insert_time << "\n";

    return new std::tuple<std::vector<size_t>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>>(min_count_list, recalls, qps, search_not_deleted, search_deleted);
}

vector<std::priority_queue<std::pair<float, labeltype>>> load_calculated_gt(vector<std::priority_queue<std::pair<float, labeltype>>> answers, std::string loc_calc_gt, size_t qsize)
{
    (vector<std::priority_queue<std::pair<float, labeltype>>>(qsize)).swap(answers);

    std::cout << "getting gt from " << loc_calc_gt << "\n";
    ifstream inputgt(loc_calc_gt);
    std::string temp_string;
    std::string dist;
    std::string labels;
    int load_qsize;
    int load_pqsize;

    inputgt >> temp_string;
    inputgt >> load_qsize;
    if (load_qsize == qsize && temp_string == "start")
    {
        inputgt >> load_pqsize;

        for (int i = 0; i < qsize; i++)
        {
            inputgt >> dist;
            inputgt >> labels;
            vector<std::string> dist_del;
            vector<std::string> labels_del;
            string temp;
            char delim = ' ';
            stringstream dd(dist);
            stringstream dl(labels);
            while (getline(dd, temp, delim))
            {
                dist_del.push_back(temp);
            };
            while (getline(dl, temp, delim))
            {
                labels_del.push_back(temp);
            };
            for (int j = 0; j < labels_del.size(); j++)
            {
                // std::cout << "check stoi\n";
                answers[i].push(std::pair<float, labeltype>(stof(dist_del[j]), (labeltype)stoi(labels_del[j])));
            }
        }
    }
    else
    {
        std::cout << "Invalid file\n";
        exit(0);
    }

    //   for(int i = 0 ; i < qsize; i++){
    //    std::cout << i << ": " << answers[i].top().first << " and " << answers[i].top().second << "\n";
    // }

    return answers;
}

void save_calculated_gt(vector<std::priority_queue<std::pair<float, labeltype>>> &answers, std::string loc_calc_gt, size_t qsize)
{
    vector<std::priority_queue<std::pair<float, labeltype>>> temp_answers(qsize);

    temp_answers.swap(answers);
    std::cout << "outputting gt to " << loc_calc_gt << "\n";
    ofstream outputgt(loc_calc_gt, ios::out | ios::trunc);
    outputgt << "start" << endl;
    outputgt << temp_answers.size() << endl;
    outputgt << temp_answers[0].size() << endl;

    for (int i = 0; i < qsize; i++)
    {
        std::priority_queue<std::pair<float, labeltype>> temp_queue;
        std::string dist = "";
        std::string labels = "";
        while (!temp_answers[i].empty())
        {
            std::pair<float, labeltype> temp_pair = temp_answers[i].top();
            temp_answers[i].pop();
            dist = dist + std::to_string(temp_pair.first) + " ";
            labels = labels + std::to_string(temp_pair.second) + " ";
            // std::cout << "pushing to answers at " << i << "\n";
            answers[i].push(temp_pair);
        }
        // std::cout << "dist: " << dist << "\n";
        // std::cout << "labels: " << labels << "\n";
        // outputgt << answers. << endl;stest_clloec
        outputgt << dist << endl;
        outputgt << labels << endl;
    }
    outputgt << "end";
    outputgt.close();
}

static void test_collection(unsigned int *massQA, unsigned char *massQ, size_t vecsize, size_t qsize, SpaceInterface<float> *l2space,
                            size_t vecdim, bool similarity, PQbinDoubleGraph<float> *appr_alg, bool use_calc_gt = false, bool save_calc_gt = false, bool load_calc_gt = false, std::string loc_calc_gt = "", bool knn_cutoff_applied = false, float knn_cutoff = 1.0, std::string index_file = "", bool hamming = false, bool weighted = false, std::string query_type = "LIGS")
{
    // std::cout << "Is labe1 1 a label? " << appr_alg->isLabel(1) << "\n";
    // std::cout << "Is labe1 899999 a label? " << appr_alg->isLabel(899999) << "\n";
    // std::cout << "Is labe1 944444 a label? " << appr_alg->isLabel(944444) << "\n";std::cout << "Begin testing\n";
    std::cout << "Build time: " << appr_alg->build_time << " s\n";
    std::cout << "Delete time: " << appr_alg->delete_time << " s\n";
    std::cout << "load calc gt: " << load_calc_gt;
    vector<std::priority_queue<std::pair<float, labeltype>>> answers;
    int k[5] = {1, 5, 10, 20, 100};
    for (int i = 0; i < 5; i++)
    {
        std::cout << "k: " << k[i] << "\n";
        if (use_calc_gt)
        {
            if (load_calc_gt)
            {
                answers = load_calculated_gt(answers, loc_calc_gt, qsize);
            }
            else
            {
                std::cout << "Calculating gt:\n";
                calc_gt(massQ, vecsize, qsize, appr_alg, l2space, vecdim, answers);
                if (save_calc_gt)
                {
                    save_calculated_gt(answers, loc_calc_gt, qsize);
                }
            }
        }
        else
        {
            std::cout << "Parsing gt:\n";
            get_gt(massQA, massQ, vecsize, qsize, l2space, vecdim, answers, k[i]);
        }
        test_vs_recall(massQ, vecsize, qsize, *appr_alg, vecdim, answers, k[i], similarity, knn_cutoff_applied, knn_cutoff, index_file, hamming, weighted, query_type);
        std::cout << "Actual memory usage: " << getCurrentRSS() / vecsize << " Mb \n";
        if (!similarity && i == 0)
        {
            break;
        }
    }
}

inline bool exists_test(const std::string &name)
{
    ifstream f(name.c_str());
    return f.good();
}

inline bool exists_test_di(const std::string &name, const int iteration_number)
{
    std::string path = "/" + name;
    for (int i = 0; i < iteration_number; i++)
    {
        std::cout << "delete_list files not present - " << path << "/" << i << "\n";
        // if(!std::filesystem::exists(path + "/" + to_string(i))){return false;} /* */
        if (!std::filesystem::exists(name + "/" + to_string(i)))
        {
            return false;
        } /* */
    }
    return true; // TODO check if working
}

std::pair<float, float> average_stdev(std::vector<float> const &v)
{
    if (v.empty())
    {
        return pair<float, float>(0, 0);
    }

    auto const count = static_cast<float>(v.size());
    float average = std::reduce(v.begin(), v.end()) / count;

    std::vector<float> diff(count);
    std::transform(v.begin(), v.end(), diff.begin(), [average](double x)
                   { return x - average; });
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / v.size());

    return pair<float, float>(average, stdev);
}

unsigned int *load_gt(std::string path_gt, int qsize)
{
    ifstream inputGT(path_gt, ios::binary);
    unsigned int *massQA = new unsigned int[qsize * 100 * 4];
    std::cout << "Loading GT...\n";
    for (int i = 0; i < qsize; i++)
    {
        int t;
        inputGT.read((char *)&t, 4);
        // std::cout << t << "\n";

        inputGT.read((char *)(massQA + 100 * i * 4), t * 4);
        /*
        for(int j = 0;j < t; j++){
            inputGT.read((char *) (massQA + 100 * i * 4 + j * 4), 4);
        }*/
        // std::cout << (char *) (massQA + 100 * i * 4) << "\n";
        if (t != 100)
        {
            std::cout << "err";
            return nullptr;
        }
    }
    inputGT.close();
    return massQA;
};

unsigned char *load_queries(std::string path_q, int qsize, int vecdim)
{
    unsigned char *massQ = new unsigned char[qsize * vecdim * 4];
    unsigned char *massb = new unsigned char[vecdim * 4];
    std::cout << "Loading queries...\n";
    ifstream inputQ(path_q, ios::binary);

    for (int i = 0; i < qsize; i++)
    {
        int in = 0;
        inputQ.read((char *)&in, 4);
        // std::cout << in << "\n";
        if (in != vecdim)
        {
            std::cout << "file error";
            exit(1);
        }

        inputQ.read((char *)massb, in * 4);
        for (int j = 0; j < vecdim * 4; j++)
        {
            massQ[i * vecdim * 4 + j] = massb[j];
        }
    }
    inputQ.close();
    std::cout << "finished loading queries\n";
    return massQ;
};

float build_standard(PQbinDoubleGraph<float> *appr_alg, std::string path_data, bool cutoff_applied, float cutoff, int vecsize, int vecdim, int mode, int high_in = -1, int low_in = 0, std::set<size_t> *nodes_to_delete = nullptr, int insert_label = 0, bool bin = false, bool calc_dist = false, std::string dist_mult = "euclidean", std::string learn_data = "", bool rand_centroids = true, bool initialize_graphs = false, float graph_cutoff = 0.5, std::vector<unsigned char *> *preloaded_data = nullptr)
{
    std::vector<float> build_times(appr_alg->max_elements_, -1);
    bool preloaded = false;
    if (path_data == "preload")
    {
        preloaded = true;
        std::cout << "preload\n";
    }
    std::vector<double> ntd_vector(1);
    ifstream input(path_data, ios::binary);
    int in = 0;
    std::cout << "dist_mult: " << dist_mult << "\n";
    unsigned char *massb = new unsigned char[vecdim * 4];
    int j1 = -1;
    std::vector<float> bin_values;
    int report_every = vecsize / 10.0;
    StopW stopwb = StopW();
    StopW stopwbf = StopW();
    int high;
    if (high_in == -1)
    {
        high = vecsize;
    }
    else
    {
        high = high_in;
    }
    int low = low_in;
    int int_add;
    if (mode == STANDARD)
    {
        int_add = vecsize;
    }
    else if (mode == PARTIAL || mode == ADDITIONAL || mode == ADDORIG || mode == ADDSIMPLE || mode == REPLACE || mode == REPLACENORECON || mode == REPLACEORIG || mode == REPLACEZERO)
    {
        int_add = high;
        if (mode == REPLACE || mode == REPLACENORECON || mode == REPLACEORIG || mode == REPLACEZERO || mode == REPLACEADDPOINT)
        {
            // std::cout << "in REPLACE\n";
            ntd_vector.resize((nodes_to_delete->size()));
            std::copy(nodes_to_delete->begin(), nodes_to_delete->end(), ntd_vector.begin());
        }
    }

    if (bin && appr_alg->bin_calculated == 0)
    {
        bin_values.reserve(vecdim * (vecsize));
    }

    std::cout << "Standard build active\n";
    std::cout << "int_add: " << int_add << " high: " << high << " low: " << low << " add_mode: " << mode << "\n";
    // std::cout << "enterpoint node: " << appr_alg->enterpoint_node_ << "\n";

    if (mode != INITIALIZE)
    {
        std::cout << "ini\n";
#pragma omp parallel for
        for (int i = 0; i < int_add; i++)
        {

            StopW stopindiv = StopW();
            unsigned char mass[vecdim * 4];
            int j2 = 0;
#pragma omp critical
            {
                if (preloaded)
                {
                    memcpy(mass, preloaded_data->at(j1 + 1), vecdim * 4);
                }
                else
                {
                    input.read((char *)&in, 4);
                    if (in != vecdim)
                    {
                        std::cout << "file error";
                        exit(1);
                    }
                    input.read((char *)massb, in * 4);
                    for (int j = 0; j < vecdim * 4; j++)
                    {
                        mass[j] = massb[j];
                    }
                }
                j1++;
                j2 = j1;
                if ((mode == PARTIAL || mode == STANDARD) && j1 % report_every == 0)
                {
                    std::cout << j1 / (0.01 * vecsize) << " %, "
                              << report_every / (1000.0 * 1e-6 * stopwb.getElapsedTimeMicro()) << " kips "
                              << " Mem: "
                              << getCurrentRSS() / 1000000 << " Mb \n";
                    stopwb.reset();
                }
            }

            // normalize_angular(mass, vecdim);

            if (j2 >= low && j2 < high)
            {

                tableint cur_c;
                StopW indiv_build_time = StopW();
                switch (mode)
                {
                case REPLACE:
                    // std::cout << "ndt_vector: " << ntd_vector[j2 - low] << "j2: " << j2 << "\n";
                    if (appr_alg->isMarkedDeletedExt(ntd_vector[j2 - low]))
                    {
                        appr_alg->unmarkDelete(ntd_vector[j2 - low]);
                    }
                    // std::cout << "REPLACE\n";atePointLabel((void *) (mass), ntd_vector[j2-low]);
                    appr_alg->updatePointLabel((void *)(mass), ntd_vector[j2 - low], cutoff_applied, cutoff);
                    break;
                case REPLACEZERO:
                    if (appr_alg->isMarkedDeletedExt(ntd_vector[j2 - low]))
                    {
                        appr_alg->unmarkDelete(ntd_vector[j2 - low]);
                    }
                    // std::cout << "REPLACEzero\n";
                    appr_alg->updatePointLabel((void *)(mass), ntd_vector[j2 - low], cutoff_applied, cutoff);
                    break;
                case REPLACEORIG:
                    if (appr_alg->isMarkedDeletedExt(ntd_vector[j2 - low]))
                    {
                        appr_alg->unmarkDelete(ntd_vector[j2 - low]);
                    }
                    // std::cout << "REPLACEorig\n";
                    appr_alg->updatePointLabel((void *)(mass), ntd_vector[j2 - low], cutoff_applied, cutoff);
                    break;
                case REPLACENORECON:
                    if (appr_alg->isMarkedDeletedExt(ntd_vector[j2 - low]))
                    {
                        appr_alg->unmarkDelete(ntd_vector[j2 - low]);
                    }
                    // std::cout << "REPLACErecom\n";
                    {
                        auto search = appr_alg->label_lookup_.find(ntd_vector[j2 - low]);
                        if (search != appr_alg->label_lookup_.end())
                        {
                            memcpy(appr_alg->getDataByLabelChar(ntd_vector[j2 - low]), (void *)(mass), appr_alg->data_size_);
                        }
                        else
                        {
                            // std::cout << "Point does not exist - using addPoint";
                            appr_alg->addPoint((void *)(mass), ntd_vector[j2 - low], cutoff_applied, cutoff);
                        }
                    }
                    break;
                case ADDITIONAL:
                    cur_c = appr_alg->addPoint((void *)(mass), (size_t)j2, cutoff_applied, cutoff);
                    /*for(int m = 0; m < vecdim * 4; m = m + 4){
                        char temp1[4];
                        char temp2[4];
                        for(int n = 0; n < 4;n++){
                            temp1[n] = mass[m + n];
                        }
                        std::cout << ((float *) temp1);
                        }*/
                    break;
                case ADDORIG:
                    cur_c = appr_alg->addPoint((void *)(mass), (size_t)(j2), cutoff_applied, cutoff);
                    break;
                case ADDSIMPLE:
                    // cur_c = appr_alg->addPointTransferLinks((void *)(mass), (size_t)(j2), (size_t)(j2 - insert_label));
                    break;
                case SELECT:
                    break;
                default:
                    // std::cout << "default\n";
                    cur_c = appr_alg->addPoint((void *)(mass), (size_t)j2, cutoff_applied, cutoff);

                    break;
                }

                build_times[j2] = indiv_build_time.getElapsedTimeMicro();
                if (dist_mult == "angular")
                {
                    if (ntd_vector.size() > 1)
                    {
                        appr_alg->normalize_angular((labeltype)ntd_vector[j2 - low], vecdim);
                    }
                    else
                    {
                        // std::cout << "ntd_vec size check\n";
                        appr_alg->normalize_angular(cur_c, vecdim);
                        /*char * checking = appr_alg->getDataByInternalId(cur_c);
                        for(int m = 0; m < vecdim * 4; m = m + 4){
                            char temp1[4];
                            char temp2[4];
                            for(int n = 0; n < 4;n++){
                                temp1[n] = mass[m + n];
                                temp2[n] = checking[m + n];
                            }
                            if(*((float *) temp1) == *((float *) temp2)){
                                std::cout << "yes " << m << " " << *((float *) temp1) << " " << *((float *) temp2) << "\n";
                            }else{
                                std::cout << "no  " << m << " " << *((float *) temp1) << " " <<*((float *) temp2) << "\n";

                            }}*/
                    }
                }
            }

            // appr_alg->initializePoint((void *) (mass), (size_t) j2, vecsize, 0);
            // std::cout << "initialized " << j2 << "\n";
        }
    }

    if (mode == STANDARD || mode == PARTIAL)
    {


        std::cout << "bin 0 size: " << appr_alg->bin_lookup_[0].size() << "\n";
        std::cout << "testa\n";
        StopW stopicl = StopW();
        std::cout << "testa2\n";
        appr_alg->vecdim = vecdim;
        appr_alg->initializeCentroidList(rand_centroids, dist_mult, learn_data);
        float icl_time = 1e-6 * stopicl.getElapsedTimeMicro();
        std::cout << "bin 0 size: " << appr_alg->bin_lookup_[0].size() << "\n";
        std::cout << "testc\n";
        std::cout << "initializeCentroidList time: " << icl_time << "\n";
        appr_alg->delete_time = icl_time;
        stopicl.reset();
        {
            std::cout << "cutoff not applied\n";
            appr_alg->initializeBins();
            std::cout << "initialized bins\n";
            // appr_alg->initializeBin();
        }

        // appr_alg->initializeGraph();
        if (initialize_graphs)
        {
            std::cout << "graph intiiniznig\n";
            // appr_alg->initializeGraph(true, graph_cutoff);
            appr_alg->initializeGraph();
        }
        1e-6 * stopicl.getElapsedTimeMicro();
        float ib_time = 1e-6 * stopicl.getElapsedTimeMicro();
        std::cout << "initializeBins time: " << ib_time << "\n";
    }


    appr_alg->build_time = 1e-6 * stopwbf.getElapsedTimeMicro();
    
    std::cout << "Standard build time: " << 1e-6 * stopwbf.getElapsedTimeMicro() << " s\n";
  

    return appr_alg->build_time; // returns the number of points added @ initial standard add mode.
}

float build_corr(PQbinDoubleGraph<float> *appr_alg, std::string path_data, bool cutoff_applied, float cutoff, int vecsize, int vecdim, int mode, int insert_label, int high_in = -1, int low_in = 0, std::set<size_t> *nodes_to_delete = nullptr, int point_selection = RANDOM, float sigma = 0.0, float norm = 0.0, bool bin = false, bool calc_dist = false, std::string dist_mult = "euclidean", std::vector<unsigned char *> *preloaded_data = nullptr)
{
    bool preloaded = false;
    if (path_data == "preload")
    {
        preloaded = true;
    }
    std::string prev = std::to_string((int)(sigma * 100.0 + 0.5));
    std::cout << "prev: " << prev << " sigma* 100 :" << (int)(sigma * 100.0 + 0.5) << " appr_alg size: " << appr_alg->max_elements_ << "\n";
    string end_text = std::string(3 - std::min((size_t)3, prev.length()), '0') + prev;

    std::vector<double> ntd_vector(1);
    ifstream input(path_data + "_" + end_text, ios::binary);
    std::cout << "loading file: " << path_data << "_" << end_text << "\n";
    int in = 0;
    unsigned char *massb = new unsigned char[vecdim * 4];
    int j1 = -1;
    std::vector<float> bin_values;
    int report_every = vecsize / 10.0;
    StopW stopwb = StopW();
    StopW stopwbf = StopW();
    int high = high_in;
    int low = low_in;
    int int_add = high;
    if (bin && appr_alg->bin_calculated == 0)
    {
        bin_values.reserve(vecdim * (vecsize));
    }

    ntd_vector.resize((nodes_to_delete->size()));
    std::copy(nodes_to_delete->begin(), nodes_to_delete->end(), ntd_vector.begin());

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, norm * sigma);

    std::cout << "correlated build active\n";
    std::cout << "int_add: " << int_add << " high: " << high << " low: " << low << " add_mode: " << mode << " sigma: " << sigma << "\n";
    // std::cout << "enterpoint node: " << appr_alg->enterpoint_node_ << "\n";
    std::cout << "ntd_vector size: " << ntd_vector.size() << "\n";

    std::cout << "ini\n";
#pragma omp parallel for
    for (int i = 0; i < int_add; i++)
    {
        unsigned char mass[vecdim * 4];
        int j2 = 0;
        StopW stopindiv = StopW();
#pragma omp critical
        {
            // std::cout << "i: " << i << "\n"; //DELETE
            if (preloaded)
            {
                if (j2 >= low && j2 < high)
                {
                    memcpy(mass, preloaded_data->at(j1 + 1), vecdim * 4);
                }
            }
            else
            {
                input.read((char *)&in, 4);
                if (in != vecdim)
                {
                    std::cout << "file error";
                    exit(1);
                }
                input.read((char *)massb, in * 4);
                for (int j = 0; j < vecdim * 4; j++)
                {
                    mass[j] = massb[j];
                }
            }
            j1++;
            j2 = j1;
            if (mode == STANDARD && j1 % report_every == 0)
            {
                std::cout << j1 / (0.01 * vecsize) << " %, "
                          << report_every / (1000.0 * 1e-6 * stopwb.getElapsedTimeMicro()) << " kips "
                          << " Mem: "
                          << getCurrentRSS() / 1000000 << " Mb \n";
                stopwb.reset();
            }
        }

        if (j2 >= low && j2 < high)
        {

            if (std::count(ntd_vector.begin(), ntd_vector.end(), j2))
            {
                tableint cur_c;

                switch (mode)
                {
                case REPLACE:
                    if (appr_alg->isMarkedDeletedExt((size_t)j2))
                    {
                        appr_alg->unmarkDelete((size_t)j2);
                    }
                    appr_alg->updatePointLabel((void *)(mass), (size_t)j2, cutoff_applied, cutoff);
                    break;
                case REPLACEZERO:
                    if (appr_alg->isMarkedDeletedExt((size_t)j2))
                    {
                        appr_alg->unmarkDelete((size_t)j2);
                    }
                    appr_alg->updatePointLabel((void *)(mass), (size_t)j2, cutoff_applied, cutoff);
                    break;
                case REPLACEORIG:
                    if (appr_alg->isMarkedDeletedExt((size_t)j2))
                    {
                        appr_alg->unmarkDelete((size_t)j2);
                    }
                    appr_alg->updatePointLabel((void *)(mass), (size_t)j2, cutoff_applied, cutoff);
                    break;
                case REPLACENORECON:
                    if (appr_alg->isMarkedDeletedExt((size_t)j2))
                    {
                        appr_alg->unmarkDelete((size_t)j2);
                    }
                    {
                        auto search = appr_alg->label_lookup_.find((size_t)j2);
                        if (search != appr_alg->label_lookup_.end())
                        {
                            memcpy(appr_alg->getDataByLabelChar((size_t)j2), (void *)(mass), appr_alg->data_size_);
                        }
                        else
                        {
                            appr_alg->addPoint<true>((void *)(mass), (size_t)j2, cutoff_applied, cutoff);
                        }
                    }
                    break;
                case ADDITIONAL:
                    cur_c = appr_alg->addPoint<true>((void *)(mass), (size_t)(j2 + insert_label), cutoff_applied, cutoff);
                    break;
                case ADDORIG:
                    cur_c = appr_alg->addPoint<true>((void *)(mass), (size_t)(j2 + insert_label), cutoff_applied, cutoff);
                    break;
                case ADDSIMPLE:
                    break;
                case SELECT:
                    break;
                default:
                    cur_c = appr_alg->addPoint<true>((void *)(mass), (size_t)j2, cutoff_applied, cutoff);
                    break;
                }
                if (dist_mult == "angular")
                {
                    if (mode == ADDITIONAL || mode == ADDORIG || mode == ADDSIMPLE)
                    {
                        appr_alg->normalize_angular((labeltype)(j2 + insert_label), vecdim);
                    }
                    else
                    {
                        appr_alg->normalize_angular((labeltype)(j2), vecdim);
                    }
                }
            }
        }

    }
    appr_alg->build_time = 1e-6 * stopwbf.getElapsedTimeMicro();
    std::cout << "Standard build time: " << 1e-6 * stopwbf.getElapsedTimeMicro() << " s\n";
    appr_alg->point_selection = point_selection;
    appr_alg->sigma = sigma;
    return appr_alg->build_time; // returns the number of points added @ initial standard add mode.
}

float delete_nodes(PQbinDoubleGraph<float> *appr_alg, int vecsize, int vecdim, std::set<size_t> *deleteListSet, int delete_mode, float neighbor_update_chance = 0.0, float percentage_dq = 0.00, int KNN = 10)
{
    std::vector<float> delete_times(appr_alg->max_elements_, -1);
    std::cout << "vecsize: " << vecsize << "\n";
    std::vector<size_t> *deleteList;
    deleteList = new std::vector<size_t>(deleteListSet->begin(), deleteListSet->end());
    bool batch_ = false;
    bool parallel_ = false;
    bool level_zero = false;
    bool delete_lazy = false;
    bool drop_deleted = false;
    std::vector<tableint> delete_list_intId(deleteList->size());

    if (delete_mode == MARKDEL)
    {
    }
    else if (delete_mode == MDUPDATE)
    {
        batch_ = true;
        parallel_ = true;
    }
    else if (delete_mode == MDUPPRUNE)
    {
        batch_ = true;
        parallel_ = true;
        drop_deleted = true;
    }
    else if (delete_mode == LAZY)
    {
        batch_ = true;
        parallel_ = true;
        delete_lazy = true;
    }
    else if (delete_mode == MDUPDATEND)
    {
        batch_ = true;
    }
    else if (delete_mode == MDUPPRUNEND)
    {
        batch_ = true;
        drop_deleted = true;
    }

    StopW stopw_del0 = StopW();

    std::cout << "deletelist_length: " << deleteList->size() << "\n";
    {
        StopW MDW = StopW();
        std::cout << "std deletion\n";
#pragma omp parallel for
        for (int i = 0; i < deleteList->size(); i++)
        {
            StopW indiv_delete_time = StopW();
            tableint deleteid;
            deleteid = appr_alg->deletePoint((size_t)deleteList->at(i));
            delete_list_intId[i] = deleteid;
            delete_times[deleteList->at(i)] = indiv_delete_time.getElapsedTimeMicro();
            // std::cout << "Deleted " << deleteList[i] << "\n";
            // std::cout << "What 4\n";
        }
        float MDTime = MDW.getElapsedTimeMicro() * 1e-6 / deleteList->size();
#pragma omp parallel for
        for (int i = 0; i < delete_list_intId.size(); i++)
        {
            appr_alg->delete_time_indiv[delete_list_intId[i]] = MDTime;
        }
    }
    appr_alg->delete_time = 1e-6 * stopw_del0.getElapsedTimeMicro();
    std::cout << "Delete time:" << 1e-6 * stopw_del0.getElapsedTimeMicro() << "  seconds\n";
    


    return appr_alg->delete_time;
}

std::set<size_t> *get_nodes_to_delete(PQbinDoubleGraph<float> *appr_alg, int vecsize, int vecdim, int high_in, int low_in, int delete_selection, float proportion, int correl_no = 1, float correl_chance = 1.0)
{
    std::set<size_t> *nodes_to_delete = new std::set<size_t>;
    int number_to_delete = (int)((((float)high_in) - ((float)low_in)) * proportion);
    std::cout << "high_in: " << high_in << " low_in: " << low_in << " proportionh: " << proportion << "\n";
    std::cout << "Selecting " << number_to_delete << " nodes to delete\n";
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed); // Mersenne Twister RNG

    // Generate a random number between 0 and 99
    std::uniform_int_distribution<int> distribution(0, INT_MAX);
    if (delete_selection == FIFO)
    {
        for (int i = low_in; i < (low_in + number_to_delete); i++)
        {
            nodes_to_delete->insert((size_t)i);
        }
    }
    else if (delete_selection == LIFO)
    {
        for (int i = high_in - number_to_delete; i < high_in; i++)
        {
            nodes_to_delete->insert((size_t)i);
        }
    }
    else if (delete_selection == CORREL)
    {
        while (nodes_to_delete->size() < number_to_delete)
        {
            int cent_node;
            cent_node = distribution(generator) % (high_in - low_in) + low_in;
            while (!(appr_alg->isLabel(cent_node)))
            {
                cent_node = distribution(generator) % (high_in - low_in) + low_in;
            }
            // std::cout << "Searching on " << cent_node << " number of neighbors: " << correl_no << "\n";
            std::priority_queue<std::pair<float, hnswlib::labeltype>> nn_list = appr_alg->searchKnn(appr_alg->getDataByLabelChar(cent_node), correl_no);

            while (!nn_list.empty() && nodes_to_delete->size() < number_to_delete)
            {
                nodes_to_delete->insert(nn_list.top().second);
                nn_list.pop();
            }
            // std::cout << nodes_to_delete->size() << std::endl;
        }
    }
    else if (delete_selection == CORRELCH)
    {
        // std::cout << "using correlch " << "get numbers to delete: " << number_to_delete << "\n";
        while (nodes_to_delete->size() < number_to_delete)
        {
            int cent_node;
            cent_node = distribution(generator) % (high_in - low_in) + low_in;
            while (!(appr_alg->isLabel(cent_node)))
            {
                cent_node = distribution(generator) % (high_in - low_in) + low_in;
            }
            std::priority_queue<std::pair<float, hnswlib::labeltype>> nn_list = appr_alg->searchKnn(appr_alg->getDataByLabelChar(cent_node), correl_no);
            while (!nn_list.empty())
            {
                if (((double)distribution(generator) / (RAND_MAX)) < correl_chance)
                {
                    nodes_to_delete->insert(nn_list.top().second);
                    // std::cout << "Cent. node: " << cent_node << " node_inserted: " <<nn_list.top().second << "\n";
                }
                nn_list.pop();
            }
        }
    }
    else if (delete_selection == RANDOM)
    {
        while (nodes_to_delete->size() < number_to_delete)
        {
            int temp = distribution(generator) % (high_in - low_in) + low_in;
            if (appr_alg->isLabel(temp))
            {
                nodes_to_delete->insert(temp);
            }
        }
    }
    return nodes_to_delete;
}

std::set<size_t> *load_nodes_to_delete(std::string delete_list, int del_sel, int vecsize)
{
    std::set<size_t> *nodes_to_delete = new std::set<size_t>;
    std::cout << "Loading delete list from " << delete_list << "\n";
    std::cout << "Del_sel: " << del_sel << ", vecsize: " << vecsize << "\n";
    ifstream inputNOD(delete_list);
    std::string temp;
    size_t temp_node = -1;
    int del_sel_read;
    int vecsize_read;
    inputNOD >> del_sel_read;
    inputNOD >> vecsize_read;
    inputNOD >> temp;
    if (del_sel != del_sel_read || vecsize != vecsize_read)
    {
        std::cout << "NoD files settings do not match\n";
        exit(0);
    }
    while (true)
    {
        inputNOD >> temp_node;
        if (temp_node == -1)
        {
            break;
        }
        nodes_to_delete->insert(temp_node);
    }
    inputNOD.close();
    return nodes_to_delete;
}

void save_nodes_to_delete(std::string delete_list, std::set<size_t> *nodes_to_delete, int del_sel, int vecsize)
{
    set<size_t>::iterator itr;
    std::cout << "outputting to " << delete_list << "\n";
    ofstream outputNOD(delete_list, ios::out | ios::trunc);
    outputNOD << del_sel << endl
              << vecsize << endl
              << "check" << endl;
    for (itr = nodes_to_delete->begin(); itr != nodes_to_delete->end(); itr++)
    {
        outputNOD << *itr << endl;
    }
    outputNOD << -1;
    outputNOD.close();
}

void preload_data(std::vector<unsigned char *> *preloaded_data, std::string path_data, int vecdim, int vecsize)
{
    int j1 = -1;
    ifstream input(path_data, ios::binary);
    int in = 0;
    unsigned char *massb = new unsigned char[vecdim * 4];

    for (int i = 0; i < vecsize; i++)
    {
        unsigned char *mass = new unsigned char[vecdim * 4];
        int j2 = 0;
#pragma omp critical
        {
            input.read((char *)&in, 4);
            if (in != vecdim)
            {
                std::cout << "file error";
                delete[] mass; // Clean up the allocated memory before exiting
                delete[] massb;
                exit(1);
            }

            input.read((char *)massb, in * 4);
            for (int j = 0; j < vecdim * 4; j++)
            {
                mass[j] = massb[j];
            }

            j1++;
            j2 = j1;

            // Insert the mass pointer to the vector
            preloaded_data->push_back((unsigned char *)mass);
        }
    }
    delete[] massb; // Release the memory after the loop
}

void delete_insert_nodes(PQbinDoubleGraph<float> *appr_mod, SpaceInterface<float> *distspace, bool cutoff_applied, float cutoff, std::string path_data, std::string path_q, std::string path_gt, int vecsize, int vecdim, float proportion, int iteration_number, bool ldl, bool sdl, std::string delete_list, int delete_mode, int update_mode, int del_sel, int del_k, float chance, std::string index_modified, bool overwrite, bool calc_gt_true, bool save_calc_gt, bool load_calc_gt, std::string loc_calc_gt, int query_steps, int qsize, int point_selection, float sigma, float norm, bool bin = false, bool calc_dist = false, std::string dist_metric = "euclidean", bool knn_cutoff_applied = false, float knn_cutoff = 1.0, bool weighted = false, std::string query_type = "LIGS")
{
    bool current_label = true;
    bool prev_label = appr_mod->isLabel(((int)vecsize * (1 - proportion) * 0.9) - 1);
    size_t insert_label;
    int num_per_iter = (int)(vecsize * proportion / iteration_number);
    int add_mode = update_mode;
    set<int> checkpoint;
    std::cout << "no of elements: " << appr_mod->traverseThroughTables() << "\n";
    if (iteration_number < 200)
    {
        checkpoint.insert(10);
        checkpoint.insert(20);
        checkpoint.insert(30);
        checkpoint.insert(40);
        checkpoint.insert(50);
        checkpoint.insert(60);
        checkpoint.insert(70);
        checkpoint.insert(80);
        checkpoint.insert(90);
    }
    checkpoint.insert(100);
    checkpoint.insert(900);
    checkpoint.insert(1000);

    std::vector<unsigned char *> *preloaded_data = new std::vector<unsigned char *>();
    preload_data(preloaded_data, path_data, vecdim, vecsize);

    std::cout << "Checking final inserted point...\n";

    for (int i = (int)vecsize * (1 - proportion) * 0.9; i < vecsize * (1 - proportion) * 1.1; i++)
    {
        current_label = appr_mod->isLabel(i);
        if (!current_label && !prev_label)
        {
            insert_label = i - 1;
            break;
        }
        prev_label = current_label;
    }

    std::cout << "Final inserted point: " << insert_label << "\n";


    float total_delete_time = 0.0;
    float total_delete_insert_time = 0.0;
    float total_insert_time = 0.0;

    int lower_point = insert_label;
    int upper_point = insert_label + num_per_iter;
    std::vector<std::string> out;
    std::cout << delete_list << "\n";
    tokenize(index_modified, '.', out);
    std::cout << "making directoy: " << delete_list << "\n";
    if (!std::filesystem::exists(delete_list))
    {
        std::filesystem::create_directory(delete_list);
    }
    std::cout << "director made\n";

    StopW itr_time = StopW();
    for (int iteration = 1; iteration <= iteration_number; iteration++)
    {
        // std::cout << iteration << "\n";
        std::streambuf *psbuf, *backup;
        if (!checkpoint.count(iteration))
        {
            backup = std::cout.rdbuf();
            std::cout.rdbuf(nullptr);
        }
        std::cout << "On " << iteration << "th loop, num_per_iter: " << num_per_iter << " upper_point: " << upper_point << " lower_point: " << lower_point << "\n";
        std::set<size_t> *nodes_to_delete;
        if (ldl)
        {
            std::cout << "Loading index from " << delete_list + "/" + to_string(iteration) << "\n";
            nodes_to_delete = load_nodes_to_delete(delete_list + "/" + to_string(iteration), del_sel, vecsize);
        }
        else
        {
            std::cout << "Making delete_node list for ";
            if (del_sel == FIFO)
            {
                std::cout << "FIFO\n";
                nodes_to_delete = get_nodes_to_delete(appr_mod, vecsize, vecdim, num_per_iter * (iteration), num_per_iter * (iteration - 1), FIFO, 1.0);
            }
            else if (del_sel == LIFO)
            {
                std::cout << "LIFO\n";
                nodes_to_delete = get_nodes_to_delete(appr_mod, vecsize, vecdim, lower_point, lower_point - num_per_iter, LIFO, 1.0);
            }
            else
            {
                std::cout << "other deletion method\n";
                if (del_sel == CORREL)
                {
                    chance = 1.0;
                }
                std::cout << lower_point << ", " << del_sel << ", " << (float)num_per_iter / (float)lower_point << ", " << del_k << ", " << chance << std::endl;
                nodes_to_delete = get_nodes_to_delete(appr_mod, vecsize, vecdim, lower_point, 0, del_sel, (float)num_per_iter / (float)lower_point, del_k, chance);
            }
        }
        // std::cout << "Length of nodes_to_delete: " << nodes_to_delete->size() << " (cf. num_per_iter: " << num_per_iter << "\n";

        if (sdl)
        {
            std::cout << "Saving index for " << delete_list + "/" + to_string(iteration) << "\n";
            if (overwrite || exists_test(delete_list + "/" + to_string(iteration)))
            {
                save_nodes_to_delete(delete_list + "/" + to_string(iteration), nodes_to_delete, del_sel, vecsize);
            }
            else
            {
                std::cout << "delete node file not saved: file exists\n";
            }
        }

        std::cout << "Deleting nodes\n";
        total_delete_time = total_delete_time + delete_nodes(appr_mod, vecsize, vecdim, nodes_to_delete, delete_mode);

        // totalElements = 0;
        // for (const auto& map : appr_mod->bin_lookup_) {
        //     for (const auto& pair : map) {
        //         totalElements += std::distance(pair.second.begin(), pair.second.end());
        //     }
        // }

        // std::cout << "Total number of elements in all forward_lists: " << totalElements << std::endl;

        std::cout << "add_mode: " << add_mode << "\n";
        if (point_selection == RANDOM)
        {

            std::cout << "Adding nodes from " << lower_point << " to " << upper_point << " in mode " << add_mode << "\n";
            std::cout << "In point selection: random\n";
            total_insert_time = total_insert_time + build_standard(appr_mod, "preload", cutoff_applied, cutoff, vecsize, vecdim, ADDITIONAL, upper_point, lower_point, nodes_to_delete, insert_label, bin, calc_dist, dist_metric, "", true, false, 0.5, preloaded_data);
        }
        else if (point_selection == CORREL)
        {

            std::cout << "Adding nodes from " << lower_point - insert_label << " to " << upper_point - insert_label << " in mode " << add_mode << "\n";
            std::cout << "In point selection: correl\n";
            total_insert_time = total_insert_time + build_corr(appr_mod, "preload", cutoff_applied, cutoff, vecsize, vecdim, ADDITIONAL, insert_label, upper_point - insert_label, lower_point - insert_label, nodes_to_delete, point_selection, sigma, norm, bin, calc_dist, dist_metric, preloaded_data);
        }
        else
        {
            std::cout << "Invaild point selection.\n";
            exit(0);
        }

        total_delete_insert_time = total_delete_time + total_insert_time;

        std::string temp;

        std::cout << "query_steps: " << query_steps << " iteration: " << iteration << " in checkpoint? " << (checkpoint.find(iteration) != checkpoint.end()) << "\n";
        // std::cin >> temp;

        vector<std::priority_queue<std::pair<float, labeltype>>> answers;
        if (checkpoint.find(iteration) != checkpoint.end())
        {
            std::cout << "process time:" << 1e-6 * itr_time.getElapsedTimeMicro() << "  seconds\n";
            appr_mod->build_time = total_insert_time;
            appr_mod->delete_time = total_delete_time;
            std::cout << "Build time: " << appr_mod->build_time << " s\n";
            std::cout << "Delete time: " << appr_mod->delete_time << " s\n";
        }
        if (query_steps == ALL || (query_steps == CHECKPOINT && checkpoint.find(iteration) != checkpoint.end()))
        {
            // std::cout << "Inside\n";

            std::cout << "total_delete_insert_time:" << total_delete_insert_time << "\n";
            unsigned char *massQ = new unsigned char[qsize * vecdim * 4];
            massQ = load_queries(path_q, qsize, vecdim);
            if (calc_gt_true)
            {

                std::cout << "calc_gt_true\n";
                if ((load_calc_gt || save_calc_gt) && !std::filesystem::exists(loc_calc_gt))
                {
                    std::filesystem::create_directory(loc_calc_gt);
                }

                if (load_calc_gt)
                {
                    answers = load_calculated_gt(answers, loc_calc_gt + "/" + to_string(iteration), qsize);
                }
                else
                {
                    std::cout << "Calculating gt:\n";
                    calc_gt(massQ, vecsize, qsize, appr_mod, distspace, vecdim, answers);
                    if (save_calc_gt)
                    {
                        if (overwrite || exists_test(loc_calc_gt + "/" + to_string(iteration)))
                        {
                            save_calculated_gt(answers, loc_calc_gt + "/" + to_string(iteration), qsize);
                        }
                        else
                        {
                            std::cout << "gt file exists: not saved\n";
                        }
                    }
                }
                std::cout << "start testing vs recall " << answers.size() << "\n";
                test_vs_recall(massQ, vecsize, qsize, *appr_mod, vecdim, answers, 1, false, knn_cutoff_applied, knn_cutoff, "", false, weighted, query_type);
                std::cout << "finished testing vs recall " << answers.size() << "\n";
            }
            else
            {
                unsigned int *massQA = new unsigned int[qsize * 100 * 4];
                massQA = load_gt(path_gt, qsize);
                get_gt(massQA, massQ, vecsize, qsize, distspace, vecdim, answers, 1);
                // test_vs_recall(massQ, vecsize, qsize, *appr_mod, vecdim, answers, 1, false, knn_cutoff_applied, knn_cutoff, "", false, "", false, weighted, query_type);
            }
        }

        if (checkpoint.find(iteration) != checkpoint.end())
        {
            // auto my_cstr = "di";
            std::cout << "Saving index\n";
            appr_mod->build_time = total_insert_time;
            appr_mod->delete_time = total_delete_time;
            appr_mod->delete_method = delete_mode;
            appr_mod->delete_selection = del_sel;
            appr_mod->delete_k = del_k;
            appr_mod->delete_chance = chance;
            // appr_mod->mode = std::hash(string("di"));
            appr_mod->update_method = update_mode;
            appr_mod->iteration = iteration;
            appr_mod->insert_steps = iteration_number;
            appr_mod->point_selection = point_selection;
            appr_mod->sigma = sigma;
            std::string itr = to_string(iteration);
            std::string itr_tot = to_string(iteration_number);
            itr.insert(0, 5 - itr.length(), '0');
            itr_tot.insert(0, 5 - itr_tot.length(), '0');
            if (out.size() > 2)
            {
                itr = out[0] + "_in_" + itr + "_of_" + itr_tot;
            }
            else
            {
                itr = out[0] + "_in_" + itr + "_of_" + itr_tot + "." + out[1];
            }
            std::cout << "no of elements: " << appr_mod->traverseThroughTables() << "\n";
            std::cout << "Saving index at " << itr << "\n";
            if (!overwrite && exists_test(itr))
            {
                std::cout << "index exists, aborting...\n";
                exit(0);
            }
            appr_mod->saveIndex(itr);

            std::ofstream output_add(itr + "_add", std::ios::trunc);
            output_add << appr_mod->subdivision;

            output_add.close();
            itr_time.reset();
        }

        // DOTO: integrate total build time and delete only time to the indices

        lower_point = upper_point;
        upper_point = upper_point + num_per_iter;
        if (upper_point > vecsize)
        {
            upper_point = vecsize;
        }
        nodes_to_delete->clear();
        nodes_to_delete->~set();
        if (!checkpoint.count(iteration))
        {
            std::cout.rdbuf(backup);
        }
        std::cout << "iteration: " << iteration << "\n";
        int asd = 0;
        // size_t totalElements = 0;

        // for (const auto& map : appr_mod->bin_lookup_) {
        //     for (const auto& pair : map) {
        //         totalElements += std::distance(pair.second.begin(), pair.second.end());
        //     }
        // }

        // std::cout << "Total number of elements in all forward_lists: " << totalElements << std::endl;
    }
};

PQbinDoubleGraph<float> *build_graph(SpaceInterface<float> *distspace, SpaceInterface<float> *distspace_sub, std::vector<std::pair<std::string, std::string>> *params, std::string path_data, std::string path_q, std::string path_learn, std::string path_gt, int vecsize, int vecdim, std::string dist_metric, std::string path_index, int qsize,
                                     float norm = 0.0, bool modify = false, PQbinDoubleGraph<float> *appr_mod = nullptr, std::string index_modified = "")
{
    PQbinDoubleGraph<float> *appr_alg;

    std::string mode_ = params->at(0).second;
    float proportion = 0.0;
    int centers = 0;
    float dist_mult = 0.0;
    float delete_add = 0.0;
    int delete_size = 0;
    int division_count = 1000;
    int number_group = 50;
    bool parallel_ = false;
    bool batch_ = false;
    bool level_zero = false;
    bool delete_lazy = false;
    float per_dq = 0.0;
    int knn_count = 10;
    bool knn_cutoff_applied = false;
    float knn_cutoff = 1.0;
    bool random_ = false;
    bool double_centers = false;
    float update_chance = 0.3;
    float cutoff_multiplier = 1.0;
    bool connections = 1.0;
    int subdivision = 0;
    bool border_test = false;
    bool drop_deleted = false;
    int delete_mode = FIFO;
    int high_in;
    int hyperplane_no = 0;
    int low_in;
    int del_sel;
    int pq_table_no = 0;
    int hyperplane_table_no = 0;
    float chance = 1.0;
    int del_k = 0;
    bool ldl = false;
    bool sdl = false;
    bool notest = false;
    bool overwrite = false;
    bool calc_dist = false;
    bool cutoff_applied = false;
    bool centroids_random = true;
    float cutoff = 1.0;
    bool bin = false;
    bool build_graph_bool = false;
    bool graph_cutoff = 0.5;
    bool hamming = false;
    bool weighted = false;
    std::string query_type = "LIGS";
    float neighbor_update_chance = 0.0;
    float percentage_dq = 0.0;
    int KNN = 10;
    int update_mode = STANDARD;
    int iteration_number;
    // int add_mode = ADDITIONAL;
    bool save_calc_gt = false;
    bool load_calc_gt = false;
    bool calc_gt_true = false;
    float sigma = 0.0;
    int point_selection = RANDOM;
    int query_steps = 1000;
    std::string loc_calc_gt = "";
    std::string delete_list = "";
    std::cout << "mode: " << mode_ << "\n";
    std::cout << "params size " << params->size() << "\n";
    for (int i = 0; i < params->size(); i++)
    {
        std::string str = params->at(i).first;
        // std::cout << str << " " << str.compare("-ni") << " check at bt\n";
        if (str.compare("-proportion") == 0)
        {
            proportion = stof(params->at(i).second);
            std::cout << "proportion = " << proportion << "\n";
        }
        else if (str.compare("--point_selection") == 0 || str.compare("-ps") == 0)
        {
            std::string temp = params->at(i).second;
            if (temp.compare("random") == 0 || temp.compare("rand") == 0)
            {
                point_selection = RANDOM;
                std::cout << "Using random point selection\n";
            }
            else if (temp.compare("corr") == 0 || temp.compare("correlated") == 0)
            {
                point_selection = CORREL;
                std::cout << "Using correlated point selection\n";
            }
            else
            {
                std::cout << "Invalid point selection: using random\n";
            }
        }
        else if (str.compare("-sigma") == 0)
        {
            sigma = stof(params->at(i).second);
        }
        else if (str.compare("--table_no_pq") == 0 || str.compare("-tnp") == 0)
        {
            pq_table_no = stoi(params->at(i).second);
            std::cout << "pq_table_no applied: " << pq_table_no << "\n";
        }
        else if (str.compare("--table_no_hyperplane") == 0 || str.compare("-tnh") == 0)
        {
            hyperplane_table_no = stoi(params->at(i).second);
            std::cout << "hyperplane_table_no applied: " << hyperplane_table_no << "\n";
        }
        else if (str.compare("--subdivision") == 0 || str.compare("-subdivision") == 0 || str.compare("-subdiv") == 0)
        {
            subdivision = stoi(params->at(i).second);
            std::cout << "subdivision applied: " << subdivision << "\n";
        }
        else if (str.compare("--centroid_number") == 0 || str.compare("-cn") == 0)
        {
            centers = stoi(params->at(i).second);
            std::cout << "centroid_no applied: " << centers << "\n";
        }
        else if (str.compare("-cutoff") == 0)
        {
            cutoff = stof(params->at(i).second);
            cutoff_applied = true;
            std::cout << "cutoff applied at " << cutoff << "\n";
        }
        else if (str.compare("-calc_dist") == 0 || str.compare("--calc_dist") == 0 || str.compare("-cd") == 0)
        {
            calc_dist = true;
            std::cout << "Calculating dist beetween old and new points "
                      << "\n";
            std::cout << "Note: only meaningful for di"
                      << "\n";
        }
        else if (str.compare("-bin") == 0)
        {
            bin = true;
            std::cout << "binning values\n";
        }
        else if (str.compare("--number_hyperplane") == 0 || str.compare("--hyperplane_number") == 0 || str.compare("-hyperplane") == 0 || str.compare("-hn") == 0 || str.compare("-nh") == 0)
        {
            hyperplane_no = stoi(params->at(i).second);
            std::cout << "hyperplane_no applied: " << hyperplane_no << "\n";
        }
        else if (str.compare("--add_mode") == 0 || str.compare("-am") == 0)
        {
            std::string temp = params->at(i).second;
            if (temp.compare("std") == 0)
            {
                update_mode = ADDITIONAL;
                std::cout << "Using add_mode additional standard\n";
            }
            else if (temp.compare("su") == 0)
            {
                update_mode = ADDSIMPLE;
                std::cout << "Using add_mode simple transfer\n";
            }
            else if (temp.compare("op") == 0)
            {
                update_mode = ADDORIG;
                std::cout << "Using add_mode original point search\n";
            }
            else
            {
                update_mode = ADDITIONAL;
                std::cout << "Invalid add_mode option: using add_mode standard\n";
            }
        }
        else if (str.compare("--query_iterations") == 0 || str.compare("-qi") == 0 || str.compare("--query_steps") == 0 || str.compare("-qs") == 0)
        {
            std::string temp = params->at(i).second;
            if (temp.compare("all") == 0)
            {
                query_steps = ALL;
                std::cout << "Querying at all steps\n";
            }
            else if (temp.compare("none") == 0)
            {
                query_steps = NONE;
                std::cout << "No queries\n";
            }
            else if (temp.compare("checkpoint") == 0)
            {
                query_steps = CHECKPOINT;
                std::cout << "Querying at checkpoints\n";
            }
            else
            {
                "Invalid query selection: no queries\n";
            }
        }
        else if (str.compare("--iteration_number") == 0 || str.compare("-in") == 0)
        {
            iteration_number = stoi(params->at(i).second);
            std::cout << "Using iteration_number = " << iteration_number << "\n";
        }
        else if (str.compare("-dist_mult") == 0)
        {
            dist_mult = stof(params->at(i).second);
            std::cout << "Using dist_mult = " << dist_mult << "\n";
        }
        else if (str.compare("-overwrite") == 0)
        {
            overwrite = true;
            std::cout << "Overwriting any preexisting files\n";
        }
        else if (str.compare("-delete_add") == 0)
        {
            delete_add = stof(params->at(i).second);
            std::cout << "Using delete_add = " << delete_add << "\n";
        }
        else if (str.compare("--delete_selection") == 0 || str.compare("-ds") == 0)
        {
            std::string temp = params->at(i).second;
            if (temp.compare("fifo") == 0)
            {
                del_sel = FIFO;
                std::cout << "delete_selection: FIFO\n";
            }
            else if (temp.compare("lifo") == 0)
            {
                del_sel = LIFO;
                std::cout << "delete_selection: LIFO\n";
            }
            else if (temp.compare("random") == 0)
            {
                del_sel = RANDOM;
                std::cout << "delete_selection: RANDOM\n";
            }
            else if (temp.compare("corr") == 0)
            {
                del_sel = CORREL;
                std::cout << "delete_selection: correlated\n";
            }
            else if (temp.compare("corrch") == 0)
            {
                del_sel = CORRELCH;
                std::cout << "delete_selection: correlated chance\n";
            }
            else
            {
                std::cout << "invalid delete_selection - using FIFO\n";
                del_sel = FIFO;
            }
        }
        else if (str.compare("--save_delete_list") == 0 || str.compare("-sdl") == 0)
        {
            if (mode_ != "delete" && mode_ != "delete_insert" && mode_ != "di")
            {
                std::cout << "mode not delete: sdl ignored\n";
            }
            else
            {
                delete_list = params->at(i).second;
                sdl = true;
                std::cout << "save_delete_list active\n";
                std::cout << "saving to: " << delete_list << "\n";
            }
        }
        else if (str.compare("--load_delete_list") == 0 || str.compare("-ldl") == 0)
        {
            if (mode_ != "delete" && mode_ != "delete_insert" && mode_ != "di")
            {
                std::cout << "mode not delete: ldl ignored\n";
            }
            else
            {
                delete_list = params->at(i).second;
                ldl = true;
                std::cout << "load_delete_list active\n";
                std::cout << "loading from: " << delete_list << "\n";
            }
        }
        else if (str.compare("--delete_mode") == 0 || str.compare("-dm") == 0)
        {
            std::string temp = params->at(i).second;
            if (temp.compare("lazy") == 0)
            {
                delete_mode = LAZY;
                std::cout << "delete_mode: lazy\n";
            }
            else if (temp.compare("markdel") == 0 || temp.compare("md") == 0)
            {
                delete_mode = MARKDEL;
                std::cout << "delete_mode: markDelete\n";
            }
            else if (temp.compare("mdupdate") == 0 || temp.compare("mdu") == 0)
            {
                delete_mode = MDUPDATE;
                std::cout << "delete_mode: markDelete + update\n";
            }
            else if (temp.compare("mduprune") == 0 || temp.compare("mdup") == 0)
            {
                delete_mode = MDUPPRUNE;
                std::cout << "delete_mode: markDelete + update + prune\n";
            }
            else if (temp.compare("mdupdatend") == 0 || temp.compare("mdund") == 0)
            {
                delete_mode = MDUPDATEND;
                std::cout << "delete_mode: markDelete + update + nodup\n";
            }
            else if (temp.compare("mduprunend") == 0 || temp.compare("mdupnd") == 0)
            {
                delete_mode = MDUPPRUNEND;
                std::cout << "delete_mode: markDelete + update + prune + nodup\n";
            }
            else if (temp.compare("rebuild") == 0)
            {
                delete_mode = REBUILD;
                std::cout << "delete_mode: REBUILD\n";
            }
            else
            {
                std::cout << "invalid delete_mode - using markDelete\n";
                delete_mode = MARKDEL;
            }
        }
        else if (str.compare("--update_mode") == 0 || str.compare("-um") == 0)
        {
            std::string temp = params->at(i).second;
            if (temp.compare("standard") == 0 || temp.compare("std") == 0)
            {
                update_mode = REPLACE;
                std::cout << "update_mode = standard replace\n";
            }
            else if (temp.compare("orig_point") == 0 || temp.compare("op") == 0)
            {
                update_mode = REPLACEORIG;
                std::cout << "update_mode: orig_point\n";
            }
            else if (temp.compare("simple_update") == 0 || temp.compare("su") == 0)
            {
                update_mode = REPLACENORECON;
                std::cout << "update_mode: no update\n";
            }
            else if (temp.compare("standard_zero") == 0 || temp.compare("std0") == 0)
            {
                update_mode = REPLACEZERO;
                std::cout << "update_mode: update_no neighbor\n";
            }
            else if (temp.compare("addpoint") == 0 || temp.compare("ap0") == 0)
            {
                update_mode = REPLACEADDPOINT;
                std::cout << "update_mode: addpoint\n";
            }
            else
            {
                update_mode = REPLACE;
                std::cout << "Invaild insert mode - using standard replace\n";
            }
            std::cout << "Using division_count = " << division_count << "\n";
        }
        else if (str.compare("-division") == 0)
        {
            division_count = stoi(params->at(i).second);
            std::cout << "Using division_count = " << division_count << "\n";
        }
        else if (str.compare("-delete_chance") == 0)
        {
            chance = stof(params->at(i).second);
            std::cout << "Using correlated delete chance = " << chance << "\n";
        }
        else if (str.compare("-delete_k") == 0)
        {
            del_k = stoi(params->at(i).second);
            std::cout << "Using correlated K = " << del_k << "\n";
        }
        else if (str.compare("-number_group") == 0)
        {
            number_group = stoi(params->at(i).second);
            std::cout << "Using number_group = " << number_group << "\n";
        }
        else if (str.compare("-level_zero") == 0)
        {
            std::cout << "only checking level_zero"
                      << "\n";
            level_zero = true;
        }
        else if (str.compare("-parallel") == 0)
        {
            std::cout << "parallel mode"
                      << "\n";
            parallel_ = true;
        }
        else if (str.compare("-batch") == 0)
        {
            std::cout << "batch mode"
                      << "\n";
            batch_ = true;
        }
        else if (str.compare("--search_mode") == 0 || str.compare("-sm") == 0)
        {
            if (params->at(i).second == "hamming")
            {
                std::cout << "hamming bin mode active \n";
                hamming = true;
                query_type = "hamming";
            }
            else if (params->at(i).second == "hamming_weighted")
            {
                std::cout << "hamming weightedbin mode active \n";
                hamming = true;
                weighted = true;
                query_type = "hamming_weighted";
            }
            else if (params->at(i).second == "LIGS_weighted")
            {
                std::cout << "hamming weighted SearchBaseLayerST  mode active \n";
                weighted = true;
                query_type = "LIGS_weighted";
            }
            else if (params->at(i).second == "LIGS_cutoff")
            {
                std::cout << "LIGS_cutoff  mode active \n";
                query_type = "LIGS_cutoff";
            }
            else if (params->at(i).second == "graph")
            {
                std::cout << "graph  mode active \n";
                query_type = "graph";
            }
            else
            {
                std::cout << "invalid search mode: SearchBaseLayerST used\n";
            }
        }
        else if (str.compare("--no_test") == 0 || str.compare("-nt") == 0)
        {
            std::cout << "no test"
                      << "\n";
            notest = true;
        }
        else if (str.compare("--percentage_dummy_queries") == 0 || str.compare("-p_dq") == 0)
        {
            per_dq = stof(params->at(i).second);
            std::cout << "percentage_dummy_queries = " << per_dq << "\n";
        }
        else if (str.compare("-KNN_count") == 0)
        {
            knn_count = stoi(params->at(i).second);
            std::cout << "Using knn_count = " << knn_count << "\n";
        }
        else if (str.compare("--update_chance") == 0 || str.compare("-uc") == 0)
        {
            update_chance = stof(params->at(i).second);
            std::cout << "Using update_chance  = " << update_chance << "\n";
        }
        else if (str.compare("-random") == 0)
        {
            random_ = true;
            std::cout << "Using randomized partitions\n";
        }
        else if (str.compare("--double_centers") == 0 || str.compare("-dkm") == 0)
        {
            double_centers = true;
            std::cout << "Using double k means\n";
        }
        else if (str.compare("-cutoff") == 0)
        {
            cutoff_multiplier = stof(params->at(i).second);
            std::cout << "Using cutoff_multiplier  = " << cutoff_multiplier << "\n";
        }
        else if (str.compare("-connections") == 0)
        {
            connections = true;
            std::cout << "Using connections\n";
        }
        else if (str.compare("-centroids") == 0)
        {
            if (params->at(i).second == "trained" || params->at(i).second == "train")
            {
                centroids_random = false;
                std::cout << "non-random centroids active\n";
            }
            else
            {
                std::cout << "random centroids\n";
            }
        }
        else if (str.compare("-border_test") == 0)
        {
            border_test = true;
            std::cout << "Using border_test\n";
        }
        else if (str.compare("-drop_deleted") == 0)
        {
            drop_deleted = true;
            std::cout << "Using drop_deleted\n";
        }
        else if (str.compare("-npc") == 0 || str.compare("--neighbor_update_chance") == 0)
        {
            chance = stof(params->at(i).second);
            std::cout << "Using neighbor_update_chance = " << neighbor_update_chance << "\n";
        }
        else if (str.compare("-build_graph") == 0)
        {
            build_graph_bool = true;
            std::cout << "building graph\n";
        }
        else if (str.compare("-graph_cutoff") == 0)
        {
            graph_cutoff = stof(params->at(i).second);
            std::cout << "graph cutoff:" << graph_cutoff << "\n";
        }
        else if (str.compare("-npc") == 0 || str.compare("--neighbor_update_chance") == 0)
        {
            chance = stof(params->at(i).second);
            std::cout << "Using neighbor_update_chance = " << neighbor_update_chance << "\n";
        }
        else if (str.compare("-pdq") == 0 || str.compare("--percentage_dummy_queries") == 0)
        {
            chance = stof(params->at(i).second);
            std::cout << "Using percentage_dummy_queries = " << percentage_dq << "\n";
        }
        else if (str.compare("-knn") == 0)
        {
            chance = stof(params->at(i).second);
            std::cout << "Using dummy query neibghbor number = " << KNN << "\n";
        }
        else if (str.compare("-calc_gt") == 0)
        {
            calc_gt_true = true;
            // std::cout << "Calculating gts\n";
        }
        else if (str.compare("--save_calc_gt") == 0 || str.compare("-scg") == 0)
        {
            // std::cout << "Saving calculated gt\n";
            save_calc_gt = true;
            calc_gt_true = true;
            loc_calc_gt = params->at(i).second;
        }
        else if (str.compare("--load_calc_gt") == 0 || str.compare("-lcg") == 0)
        {
            // std::cout << "Loading calculated gt\n";
            load_calc_gt = true;
            calc_gt_true = true;
            loc_calc_gt = params->at(i).second;
        }
    }
    if (modify)
    {
        std::cout << "loading index data\n";
        appr_mod->loadIndex(path_index, distspace, distspace_sub);
        std::ifstream output2(path_index + "_add2", std::ios::binary);
        readBinaryPOD(output2, subdivision);
        readBinaryPOD(output2, centers);
        readBinaryPOD(output2, hyperplane_no);
        readBinaryPOD(output2, hyperplane_table_no);
        int total_table_no;
        readBinaryPOD(output2, total_table_no);
        readBinaryPOD(output2, pq_table_no);
    }
    // td::cout <<"hcheck apa,s\n";
    StopW fulltime = StopW();
    // std::cout <<"hcheck apa,s\n";
    std::set<size_t> *nodes_to_delete = new std::set<size_t>;
    std::cout << "build new index\n";
    std::pair<float, float> build_info;
    std::ofstream out("log.txt", std::ios_base::app);
    if (hyperplane_table_no + pq_table_no == 0 && !modify)
    {

        throw std::runtime_error("No tables active");
    }
    if (hyperplane_table_no != 0 && hyperplane_no == 0 && !modify)
    {

        throw std::runtime_error("hyperplane tables active but hyperplane number = 0");
    }
    if (pq_table_no != 0 && (subdivision == 0 || centers == 0) && !modify)
    {

        throw std::runtime_error("PQ tables not zero but subdivision or centers zero");
    }
    std::cout << "total tables: " << pq_table_no + hyperplane_table_no << "\n";
    std::cout << "distspace_sub: " << distspace_sub->get_data_size() << " distspace: " << distspace->get_data_size() << "\n";
    appr_alg = new PQbinDoubleGraph<float>(distspace, distspace_sub, vecsize, vecdim, subdivision, centers, hyperplane_no, pq_table_no, hyperplane_table_no);
    std::cout << "new appr_alg, sudivision: " << appr_alg->subdivision << ", pq_table_no: " << appr_alg->pq_table_number << ", centers: " << appr_alg->centroid_no << "\n";
    std::cout << appr_alg->bin_lookup_[0].size();
    // appr_alg->mode= mode_;
    if (!modify)
    {
        std::cout << "Set up appr_alg\n";
        if (mode_ == "standard")
        {
            build_standard(appr_alg, path_data, cutoff_applied, cutoff, vecsize, vecdim, STANDARD, -1, 0, nullptr, 0, bin, calc_dist, dist_metric, path_learn, centroids_random, build_graph_bool, graph_cutoff);
        }
        if (mode_ == "partial")
        {
            std::cout << "partial construction on\n";
            std::cout << "dist_metric: " << dist_metric << "\n";
            build_standard(appr_alg, path_data, cutoff_applied, cutoff, vecsize, vecdim, PARTIAL, vecsize - (int)((float)vecsize * proportion), 0, nullptr, 0, bin, calc_dist, dist_metric, path_learn, centroids_random, build_graph_bool, graph_cutoff);
        }
        else if (mode_ == "delete")
        {
            high_in = vecsize;
            low_in = 0;
            build_standard(appr_alg, path_data, cutoff_applied, cutoff, vecsize, vecdim, STANDARD, high_in, low_in, nullptr, 0, bin, calc_dist, dist_metric, "", true, build_graph_bool, graph_cutoff);

            if (ldl == true)
            {
                std::cout << "loading nodes from " << delete_list << "...\n";
                nodes_to_delete = load_nodes_to_delete(delete_list, del_sel, vecsize);
            }
            else
            {
                nodes_to_delete = get_nodes_to_delete(appr_alg, vecsize, vecdim, high_in, low_in, del_sel, proportion, del_k, chance);
            }

            if (sdl == true)
            {
                std::cout << "saving nodes to " << delete_list << "...\n";
                save_nodes_to_delete(delete_list, nodes_to_delete, del_sel, vecsize);
            }
            else
            {
                delete_nodes(appr_alg, vecsize, vecdim, nodes_to_delete, delete_mode, neighbor_update_chance, percentage_dq, KNN);
                appr_alg->delete_method = delete_mode;
                appr_alg->delete_selection = del_sel;
                appr_alg->delete_k = del_k;
                appr_alg->delete_chance = chance;
            }
        }
        else if (mode_ == "delete_insert" || mode_ == "di")
        {
            // TODO

            std::cout << "currently delete_insert only functional with modify\n";

            // std::cout << "currently delete_insert only functional with modify\n";
            // exit(0);
        }

        else
        {
            std::cout << "invalid mode selection\n";
        }

        out.close();
        std::cout << "Time of full process: " << 1e-6 * fulltime.getElapsedTimeMicro() << " seconds\n";
        std::cout << "finished building\n";
        std::hash<std::string> hasher;
        auto hashed = hasher(mode_);
        appr_alg->mode = hashed;
        return appr_alg;
    }
    else
    {
        if (mode_ == "delete")
        {
            high_in = vecsize;
            low_in = 0;

            if (ldl == true)
            {
                std::cout << "loading nodes from " << delete_list << "...\n";
                nodes_to_delete = load_nodes_to_delete(delete_list, del_sel, vecsize);
            }
            else
            {
                nodes_to_delete = get_nodes_to_delete(appr_mod, vecsize, vecdim, high_in, low_in, del_sel, proportion, del_k, chance);
            }
            // note change form appralg to apprmod
            if (sdl == true)
            {
                if (!overwrite && exists_test(delete_list))
                {
                    std::cout << "Node list file exists. Aborting... \n";
                    exit(0);
                }
                std::cout << "saving nodes to " << delete_list << "...\n";
                save_nodes_to_delete(delete_list, nodes_to_delete, del_sel, vecsize);
            }

            else
            {
                delete_nodes(appr_mod, vecsize, vecdim, nodes_to_delete, delete_mode);

                // appr_mod->mode = mode_;
                appr_mod->delete_method = delete_mode;
                appr_mod->delete_selection = del_sel;
                appr_mod->delete_k = del_k;
                appr_mod->delete_chance = chance;
            }
        }
        else if (mode_ == "delete_insert" || mode_ == "di")
        {
            std::cout << "di active\n";
            if (appr_mod->cur_element_count > ceil(vecsize * (1 - proportion)))
            {
                std::cout << "No. of elements in currently loaded index: " << appr_mod->cur_element_count << " number of elements allowed: " << ceil(vecsize * (1 - proportion)) << "\n";
                std::cout << "Too many points in starting index! Aborting...\n";
                exit(0);
            }
            else
            {
                std::cout << "delete_insert_nodes has been activated\n";
                // std::cout << "point selection: " << point_selection << "\n";
                delete_insert_nodes(appr_mod, distspace, cutoff_applied, cutoff, path_data, path_q, path_gt, vecsize, vecdim, proportion, iteration_number, ldl, sdl, delete_list, delete_mode, update_mode, del_sel, del_k, chance, index_modified, overwrite, calc_gt_true, save_calc_gt, load_calc_gt, loc_calc_gt, query_steps, qsize, point_selection, sigma, norm, bin, calc_dist, dist_metric, knn_cutoff_applied, knn_cutoff, false, query_type);
            }
        }
        else
        {
            std::cout << "invalid mode selection for modify\n";
        }

        out.close();
        if (delete_mode == REBUILD)
        {
            appr_mod->delete_time = 1e-6 * fulltime.getElapsedTimeMicro();
        }
        if (delete_mode == REBUILD)
        {
            appr_alg->delete_time = 1e-6 * fulltime.getElapsedTimeMicro();
        }
        std::cout << "Time of full process: " << 1e-6 * fulltime.getElapsedTimeMicro() << " seconds\n";
        std::cout << "finished building\n";
        std::hash<std::string> hasher;
        auto hashed = hasher(mode_);
        appr_mod->mode = hashed;
        return appr_mod;
    }
    // partition/bootstrap/partition_percent/bootsrtap_percent/bootstrap_pq/bootstrap_pq_dist/bootstrap_pq_overlap/
}

void build_test(std::vector<std::pair<std::string, std::string>> *params, vector<size_t> *dimensions_, std::vector<std::string> *path_info)
{

    int efConstruction = 40;
    int M = 16;
    int subdiv = 1;
    bool build_new_index = false;
    bool similarity = false;
    bool save_index_only = false;
    bool save_index = false;
    bool knn_cutoff_applied = false;
    float knn_cutoff = 1.0;
    bool load_index = false;
    bool calc_gt_true = false;
    bool modify = false;
    bool notest = false;
    bool save = false;
    bool overwrite = false;
    bool save_calc_gt = false;
    bool load_calc_gt = false;
    bool bin = false;
    bool hamming = false;
    bool weighted = false;
    std::string query_type = "LIGS";
    std::string loc_calc_gt = "";
    std::string index_modified;
    std::string save_path = "";
    StopW stopw_full = StopW();
    std::cout << params->size() << " params\n";
    for (int i = 0; i < params->size(); i++)
    {
        std::string str = params->at(i).first;
        std::cout << str << "\n";
        if (str.compare("-ef") == 0)
        {
            efConstruction = stoi(params->at(i).second);
            std::cout << "Using efConstruction = " << efConstruction << "\n";
        }
        else if (str.compare("-overwrite") == 0)
        {
            overwrite = true;
        }
        else if (str.compare("-weighted") == 0)
        {
            std::cout << "check wegihted\n";
            weighted = true;
        }
        else if (str.compare("-save") == 0)
        {
            save = true;
            save_path = params->at(i).second;
        }
        else if (str.compare("-M") == 0)
        {
            M = stoi(params->at(i).second);
            std::cout << "Using M = " << M << "\n";
        }
        else if (str.compare("--new_index") == 0 || str.compare("-ni") == 0)
        {
            build_new_index = true;
            std::cout << "Building new index\n";
        }
        else if (str.compare("--similarity") == 0 || str.compare("-sim") == 0)
        {
            similarity = true;
            std::cout << "Using similarity\n";
        }
        else if (str.compare("--no_test") == 0 || str.compare("-nt") == 0)
        {
            std::cout << "no test"
                      << "\n";
            notest = true;
        }
        else if (str.compare("--knn_cutoff") == 0 || str.compare("-kc") == 0)
        {
            knn_cutoff_applied = true;
            knn_cutoff = stof(params->at(i).second);
            std::cout << "knn_cutoff applied: " << knn_cutoff << "\n";
        }
        else if (str.compare("--save_index_only") == 0 || str.compare("-sio") == 0)
        {
            save_index_only = true;
            std::cout << "Only saving index\n";
        }
        else if (str.compare("--save_index") == 0 || str.compare("-si") == 0)
        {
            save_index = true;
            std::cout << "Saving index\n";
        }
        else if (str.compare("--search_mode") == 0 || str.compare("-sm") == 0)
        {
            if (params->at(i).second == "hamming")
            {
                std::cout << "hamming bin mode active \n";
                hamming = true;
                query_type = "hamming";
            }
            else if (params->at(i).second == "hamming_weighted")
            {
                std::cout << "hamming weightedbin mode active \n";
                hamming = true;
                weighted = true;
                query_type = "hamming_weighted";
            }
            else if (params->at(i).second == "LIGS_weighted")
            {
                std::cout << "hamming weighted SearchBaseLayerST  mode active \n";
                weighted = true;
                query_type = "LIGS_weighted";
            }
            else if (params->at(i).second == "LIGS_cutoff")
            {
                std::cout << "LIGS_cutoff  mode active \n";
                query_type = "LIGS_cutoff";
            }
            else if (params->at(i).second == "graph")
            {
                std::cout << "graph  mode active \n";
                query_type = "graph";
            }
            else
            {
                std::cout << "invalid search mode: SearchBaseLayerST used\n";
            }
        }
        else if (str.compare("-bin") == 0)
        {
            bin = true;
            std::cout << "Using bin\n";
        }
        else if (str.compare("--load_index") == 0 || str.compare("-li") == 0)
        {
            load_index = true;
            std::cout << "Loading index: " << path_info->at(3) << "\n";
        }
        else if (str.compare("-calc_gt") == 0)
        {
            calc_gt_true = true;
            std::cout << "Calculating gts\n";
        }
        else if (str.compare("-modify") == 0)
        {
            std::cout << "modify active\n";
            modify = true;
            index_modified = params->at(i).second;
        }
        else if (str.compare("--save_calc_gt") == 0 || str.compare("-scg") == 0)
        {
            std::cout << "Saving calculated gt\n";
            save_calc_gt = true;
            loc_calc_gt = params->at(i).second;
            calc_gt_true = true;
        }
        else if (str.compare("--subdivision") == 0 || str.compare("-subdivision") == 0 || str.compare("-subdiv") == 0)
        {
            subdiv = stoi(params->at(i).second);
            std::cout << "subdivision applied: " << subdiv << "\n";
        }
        else if (str.compare("--load_calc_gt") == 0 || str.compare("-lcg") == 0)
        {
            std::cout << "Loading calculated gt\n";
            load_calc_gt = true;
            calc_gt_true = true;
            loc_calc_gt = params->at(i).second;
        }
    }
    float norm;
    size_t vecsize = dimensions_->at(0);
    size_t qsize = 10000; // dimensions_->at(1);
    size_t vecdim = dimensions_->at(2);
    // std::cout << "vecdim: " << vecdim << "\n";
    std::string path_gt = path_info->at(0);
    std::string path_q = path_info->at(1);
    std::string path_data = path_info->at(2);
    std::string path_index = path_info->at(3);
    std::string dist_metric = path_info->at(4);
    std::string path_learn = path_info->at(5);

    unsigned char *massbs = new unsigned char[4];

    unsigned int *massQA = new unsigned int[qsize * 100 * 4];
    massQA = load_gt(path_gt, qsize);

    unsigned char *massQ = new unsigned char[qsize * vecdim * 4];
    massQ = load_queries(path_q, qsize, vecdim);

    PQbinDoubleGraph<float> *appr_alg;
    SpaceInterface<float> *distspace;
    SpaceInterface<float> *distspace_sub;
    std::cout << "dist_metric: " << dist_metric << "\n";

    std::ifstream input(path_index + "_add2", std::ios::binary);
    readBinaryPOD(input, subdiv);
    std::cout << "vecdim : " << vecdim << " subdiv: " << subdiv << "\n";
    if (dist_metric == "euclidean")
    {
        std::cout << "new euclidean space\n";
        distspace = new L2Space(vecdim);
        distspace_sub = new L2Space(vecdim / subdiv);
    }
    else if (dist_metric == "angular" || dist_metric == "cosine")
    {
        std::cout << "new angular space\n";
        distspace = new AngularSpace(vecdim);              // done
        distspace_sub = new AngularSpace(vecdim / subdiv); // done
    }
    else if (dist_metric == "inner_product")
    {
        std::cout << "new inner product space\n";
        distspace = new InnerProductSpace(vecdim);
        distspace_sub = new InnerProductSpace(vecdim / subdiv);
    }

    if (!build_new_index && exists_test(path_index) && load_index)
    {
        std::cout << "Loading index from " << path_index << ":\n";
        appr_alg = new PQbinDoubleGraph<float>(distspace, distspace_sub, path_index, false);
        std::cout << "Actual memory usage: " << getCurrentRSS() / 1000000 << " Mb \n";
    }
    else
    {
        if (load_index)
        {
            std::cout << "no valid index found; building new index\n";
        }
        // ppr_alg = nullptr;
        // appr_alg = new PQbinDoubleGraph<float>(distspace, vecsize, M, efConstruction);
        std::cout << "vecdim : " << vecdim << " subdiv: " << subdiv << "\n";

        std::cout << "distspace_sub: " << distspace_sub->get_data_size() << " distspace: " << distspace->get_data_size() << "\n";
        appr_alg = build_graph(distspace, distspace_sub, params, path_data, path_q, path_learn, path_gt, vecsize, vecdim, dist_metric, path_index, qsize);

        std::cout << "Graph built\n";
    }
    if ((build_new_index || save_index_only || save_index) && !load_index && !modify)
    {
        if (!overwrite && exists_test(path_index))
        {
            std::cout << "Index file exists. Aborting...\n";
            exit(0);
        }
        else
        {
            appr_alg->saveIndex(path_index);
            std::ofstream output_add(path_index + "_add", std::ios::trunc);
            output_add << appr_alg->subdivision;

            output_add.close();
        }
        std::cout << "Index saved to: " << path_index << "\n";
    }
    else if (modify)
    {
        std::cout << "modify active build_test\n";
        appr_alg = build_graph(distspace, distspace_sub, params, path_data, path_q, path_learn, path_gt, vecsize, vecdim, dist_metric, path_index, qsize, norm, true, appr_alg, index_modified);
        if (!overwrite && exists_test(index_modified))
        {
            std::cout << "Index file exists. Aborting...\n";
            exit(0);
        }
        else
        {
            appr_alg->saveIndex(index_modified);
            std::ofstream output_add(index_modified + "_add", std::ios::trunc);
            output_add << appr_alg->subdivision;

            output_add.close();
        }
        std::cout << "Modified index saved to: " << index_modified << "\n";
    }

    float total_delete_time = 0.0;
    int deleted_c = 0;
    for (int i = 0; i < appr_alg->cur_element_count; i++)
    {
        if (appr_alg->isMarkedDeleted(i))
        {
            // std::cout << "i: " << i << " delete_time: " << appr_alg->delete_time_indiv[i] << "\n";
            total_delete_time = total_delete_time + appr_alg->delete_time_indiv[i];
            deleted_c++;
        }
    }
    std::cout << "total delete time: " << total_delete_time << " avage delete time: " << total_delete_time / deleted_c << " acutal delete time: " << appr_alg->delete_time << "\n";

    if (!save_index_only && !notest)
    {
        test_collection(massQA, massQ, vecsize, qsize, distspace, vecdim, similarity, appr_alg, calc_gt_true, save_calc_gt, load_calc_gt, loc_calc_gt, knn_cutoff_applied, knn_cutoff, path_index, hamming, weighted, query_type);
    }
}

constexpr unsigned int str2int(const char *str, int h = 0)
{
    return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}


////#####################################################Python Interface ######################333///

std::vector<unsigned char *> convert_numpy_matrix(py::array_t<float> input)
{
    // Check if the input is a two-dimensional array
    if (input.ndim() != 2)
    {
        throw std::runtime_error("Input should be a 2D array.");
    }

    // Prepare the vector that will hold pointers to the byte representations of rows
    std::vector<unsigned char *> preloaded_data;

    // Get the buffer info from the input array
    py::buffer_info buf_info = input.request();
    float *ptr = static_cast<float *>(buf_info.ptr);
    size_t num_rows = buf_info.shape[0];
    size_t num_cols = buf_info.shape[1];

    // Allocate and fill the vector
    for (size_t i = 0; i < num_rows; ++i)
    {
        // Allocate memory for this row's byte representation
        unsigned char *row_bytes = new unsigned char[num_cols * sizeof(float)];

        // For each element in the row, convert to bytes and copy
        for (size_t j = 0; j < num_cols; ++j)
        {
            std::memcpy(row_bytes + j * sizeof(float), &ptr[i * num_cols + j], sizeof(float));
        }

        // Add the byte array to the vector
        preloaded_data.push_back(row_bytes);
    }

    return preloaded_data;
}

PQbinDoubleGraph<float> *initialize_graph_python(int vecsize, int vecdim, int subdivision, int centers, int hyperplane_count, int pq_table_no, int hyperplane_no, std::string dist_metric)
{

    StopW fulltime = StopW();
    PQbinDoubleGraph<float> *appr_alg;

    SpaceInterface<float> *distspace;
    SpaceInterface<float> *distspace_sub;

    std::string loc_calc_gt = "";
    std::string delete_list = "";

    if (dist_metric == "euclidean")
    {
        std::cout << "new euclidean space\n";
        distspace = new L2Space(vecdim);
        distspace_sub = new L2Space(vecdim / subdivision);
    }
    else if (dist_metric == "angular" || dist_metric == "cosine")
    {
        std::cout << "new angular space\n";
        distspace = new AngularSpace(vecdim);                   // done
        distspace_sub = new AngularSpace(vecdim / subdivision); // done
    }
    else if (dist_metric == "inner_product")
    {
        std::cout << "new inner product space\n";
        distspace = new InnerProductSpace(vecdim);
        distspace_sub = new InnerProductSpace(vecdim / subdivision);
    }

    appr_alg = new PQbinDoubleGraph<float>(distspace, distspace_sub, vecsize, vecdim, subdivision, centers, hyperplane_count, pq_table_no, hyperplane_no);

    // td::cout <<"hcheck apa,s\n";
    // std::cout <<"hcheck apa,s\n";

    std::cout << "Time of graph generation: " << 1e-6 * fulltime.getElapsedTimeMicro() << " seconds\n";
    return appr_alg;
    // partition/bootstrap/partition_percent/bootsrtap_percent/bootstrap_pq/bootstrap_pq_dist/bootstrap_pq_overlap/
}

std::vector<float> get_value(PQbinDoubleGraph<float> *appr_alg, labeltype label){
    return appr_alg->getDataByLabel<float>(label);
}


float build_initial_python(PQbinDoubleGraph<float> *appr_alg, py::array_t<float> input_array, std::vector<int> &label, int vecsize, int vecdim, std::string dist_mult = "euclidean", std::string learn_data = "", bool rand_centroids = true)
{
    std::vector<unsigned char *> preloaded_data = convert_numpy_matrix(input_array);

    std::vector<float> build_times(appr_alg->max_elements_, -1);
    assert(label.size() == preloaded_data.size());
    std::cout << "dist_mult: " << dist_mult << "\n";
    std::vector<float> bin_values;
    int report_every = vecsize / 10.0;
    StopW stopwb = StopW();
    StopW stopwbf = StopW();
    bool cutoff_applied = false;
    float cutoff = 1.0;
    // if (bin && appr_alg->bin_calculated == 0)
    // {
    //     bin_values.reserve(vecdim * (vecsize));
    // }

    std::cout << "Standard build active\n";
    // std::cout << "enterpoint node: " << appr_alg->enterpoint_node_ << "\n";

    std::cout << "ini\n";
#pragma omp parallel for
    for (int i = 0; i < label.size(); i++)
    {

        StopW stopindiv = StopW();
        unsigned char mass[vecdim * 4];
#pragma omp critical
        {

            memcpy(mass, preloaded_data.at(i), vecdim * 4);

            if (i % report_every == 0)
            {
                std::cout << i / (0.01 * vecsize) << " %, "
                          << report_every / (1000.0 * 1e-6 * stopwb.getElapsedTimeMicro()) << " kips "
                          << " Mem: "
                          << getCurrentRSS() / 1000000 << " Mb \n";
                stopwb.reset();
            }
        }

        // normalize_angular(mass, vecdim);
        tableint cur_c;
        StopW indiv_build_time = StopW();
        cur_c = appr_alg->addPoint((void *)(mass), (size_t)label.at(i), cutoff_applied, cutoff);

        build_times[label.at(i)] = indiv_build_time.getElapsedTimeMicro();
        if (dist_mult == "angular")
        {

            appr_alg->normalize_angular(cur_c, vecdim);
        }
    }

    //std::cout << "test\n";
    std::cout << "bin 0 size: " << appr_alg->bin_lookup_[0].size() << "\n";
    //std::cout << "testaa\n";
    StopW stopicl = StopW();
    //std::cout << "testaa2\n";
    appr_alg->vecdim = vecdim;
    //std::cout << "testaa3\n";
    // appr_alg->initializeCentroidList(rand_centroids, dist_mult, learn_data);########################
    std::cout << "dist_mult: " << dist_mult << "\n";
    appr_alg->initializeCentroidList(true, dist_mult, "");
    //std::cout << "testaa4\n";
    // appr_alg->initializeCentroidDistanceStorage();
    float icl_time = 1e-6 * stopicl.getElapsedTimeMicro();
    //std::cout << "testaa5\n";
    std::cout << "bin 0 size: " << appr_alg->bin_lookup_[0].size() << "\n";
    std::cout << "initializeCentroidList time: " << icl_time << "\n";
    appr_alg->delete_time = icl_time;
    stopicl.reset();
    {
        std::cout << "cutoff not applied\n";
        appr_alg->initializeBins();
        std::cout << "initialized bins\n";
        // appr_alg->initializeBin();
    }

    appr_alg->build_time = 1e-6 * stopwbf.getElapsedTimeMicro();
    std::cout << "Time of build process: " << 1e-6 * stopwbf.getElapsedTimeMicro() << " seconds\n";
    std::cout << "finished building\n";
    return appr_alg->build_time; // returns the number of points added @ initial standard add mode.
}

float delete_nodes_python(PQbinDoubleGraph<float> *appr_alg, std::vector<int> &nodes_to_delete)
{
    std::vector<float> delete_times(appr_alg->max_elements_, -1);
    std::vector<size_t> *deleteList;
    std::vector<tableint> delete_list_intId(nodes_to_delete.size());

    StopW stopw_del0 = StopW();
    // std::cout << "batch: " << batch_ << " parallel_: " << parallel_ << " level_Zero: " << level_zero << " lazy: " << delete_lazy << " derop_deleted: " << drop_deleted << "\n";

    std::cout << "deletelist_length: " << nodes_to_delete.size() << "\n";
    {
        StopW MDW = StopW();
        std::cout << "std deletion\n";
#pragma omp parallel for
        for (int i = 0; i < nodes_to_delete.size(); i++)
        {
            StopW indiv_delete_time = StopW();
            tableint deleteid;
            deleteid = appr_alg->deletePoint((size_t)nodes_to_delete.at(i));
            delete_list_intId[i] = deleteid;
            delete_times[nodes_to_delete.at(i)] = indiv_delete_time.getElapsedTimeMicro();
            // std::cout << "Deleted " << deleteList[i] << "\n";
            // std::cout << "What 4\n";
        }
        float MDTime = MDW.getElapsedTimeMicro() * 1e-6 / nodes_to_delete.size();
#pragma omp parallel for
        for (int i = 0; i < delete_list_intId.size(); i++)
        {
            appr_alg->delete_time_indiv[delete_list_intId[i]] = MDTime;
        }
    }
    appr_alg->delete_time = appr_alg->delete_time + 1e-6 * stopw_del0.getElapsedTimeMicro();
    std::cout << "Delete time:" << 1e-6 * stopw_del0.getElapsedTimeMicro() << "  seconds\n";

    return appr_alg->delete_time;
}

float insert_nodes_python(PQbinDoubleGraph<float> *appr_alg, py::array_t<float> input_array, std::vector<int> &label, int vecdim, string dist_mult = "euclidean")
{
    std::vector<unsigned char *> preloaded_data = convert_numpy_matrix(input_array);

    std::vector<float> build_times(appr_alg->max_elements_, -1);
    assert(label.size() == preloaded_data.size());
    std::cout << "dist_mult: " << dist_mult << "\n";
    std::vector<float> bin_values;
    StopW stopwb = StopW();
    StopW stopwbf = StopW();


    bool cutoff_applied = false;
    float cutoff = 1.0;

    std::cout << "Standard build active\n";

    std::cout << "ini\n";
#pragma omp parallel for
    for (int i = 0; i < label.size(); i++)
    {

        StopW stopindiv = StopW();
        unsigned char mass[vecdim * 4];
#pragma omp critical
        {

            memcpy(mass, preloaded_data.at(i), vecdim * 4);
        }


        tableint cur_c;
        StopW indiv_build_time = StopW();
        cur_c = appr_alg->addPoint((void *)(mass), (size_t)label.at(i), cutoff_applied, cutoff);

        build_times[label.at(i)] = indiv_build_time.getElapsedTimeMicro();
        if (dist_mult == "angular")
        {

            appr_alg->normalize_angular(cur_c, vecdim);
        }
    }

    std::cout << "bin 0 size: " << appr_alg->bin_lookup_[0].size() << "\n";
    StopW stopicl = StopW();
    appr_alg->vecdim = vecdim;

    std::cout << "test1.5\n";
    for (auto &ptr : preloaded_data)
    {
        delete[] ptr;
    }
    std::cout << "test2\n";
    appr_alg->build_time = appr_alg->build_time + 1e-6 * stopwbf.getElapsedTimeMicro();
    std::cout << "Time of build process: " << appr_alg->build_time << " seconds\n";
    std::cout << "finished building\n";
    return appr_alg->build_time; // returns the number of points added @ initial standard add mode.
}

bool check_label(PQbinDoubleGraph<float> *appr_alg, size_t label) {
    return appr_alg->isLabel(label);

}

void convert_to_queries(py::array_t<float> input_array, unsigned char *massQ) {
    auto buf = input_array.unchecked<2>(); // Access array elements without bounds checking
    auto shape = input_array.shape();
    size_t qsize = shape[0]; // Number of rows
    size_t qdim = shape[1];  // Number of columns

    // Convert/cast each float to unsigned char
    for (size_t j = 0; j < qsize; ++j) {
        for (size_t i = 0; i < qdim; ++i) {
            float value = buf(j, i); // Get value as float
            unsigned char* bytes = reinterpret_cast<unsigned char*>(&value);
            // Copying the float bytes to the massQ buffer
            std::copy(bytes, bytes + sizeof(float), massQ + (j * qdim + i) * sizeof(float));
        }
    }
    // No return needed since massQ is modified directly
}

std::vector<std::vector<std::pair<float, labeltype>>> convertPriorityQueuesToVectors(
    const std::vector<std::priority_queue<std::pair<float, labeltype>>>& queues) {
    std::vector<std::vector<std::pair<float, labeltype>>> result;
    for (const auto& queue : queues) {
        std::vector<std::pair<float, labeltype>> vec;
        auto temp_queue = queue; // Copy the queue to not modify the original
        while (!temp_queue.empty()) {
            vec.push_back(temp_queue.top());
            temp_queue.pop();
        }
        std::reverse(vec.begin(), vec.end()); // Optional: reverse to maintain the priority queue's order
        result.push_back(vec);
    }
    return result;
}

std::vector<std::priority_queue<std::pair<float, labeltype>>> convertVectorsToPriorityQueues(
    const std::vector<std::vector<std::pair<float, labeltype>>>& vecs) {
    std::vector<std::priority_queue<std::pair<float, labeltype>>> queues;
    for (const auto& vec : vecs) {
        std::priority_queue<std::pair<float, labeltype>> pq;
        for (const auto& pair : vec) {
            pq.push(pair);
        }
        queues.push_back(pq);
    }
    return queues;
}

std::vector<std::priority_queue<std::pair<float, labeltype>>> calculate_gt_python(PQbinDoubleGraph<float> *appr_alg, py::array_t<float> query_array, string loc_save_gt, string dist_metric = "euclidean")
{

    auto shape = query_array.shape();
    size_t qsize = shape[0]; // Number of rows
    size_t qdim = shape[1];  // Number of columns
    std::cout << "Qsize: " << qsize << " ,qdim: " << qdim << "\n";
    size_t vecsize = appr_alg->max_elements_;
    size_t vecdim = qdim;

    // Allocate memory for massQ
    unsigned char* massQ = new unsigned char[qsize * qdim * sizeof(float)];

    // Convert the query array to massQ
    convert_to_queries(query_array, massQ);
    std::cout << "massQ converted\n";
    SpaceInterface<float> *distspace;
    if (dist_metric == "euclidean")
    {
        distspace = new L2Space(qdim);
    }
    else if (dist_metric == "angular" || dist_metric == "cosine")
    {
        distspace = new AngularSpace(qdim);
    } // done
    else if (dist_metric == "inner_product")
    {
        distspace = new InnerProductSpace(qdim);
    }
    std::vector<std::priority_queue<std::pair<float, labeltype>>> answers(qsize);
    std::cout << "to calc_gt\n";
    calc_gt(massQ, vecsize, qsize, appr_alg, distspace, qdim, answers);
    std::cout << "saving\n";
    save_calculated_gt(answers, loc_save_gt, qsize);
    return answers;
}



std::vector<std::priority_queue<std::pair<float, labeltype>>> load_gt_python(string loc_save_gt, size_t qsize)
{
    std::vector<std::priority_queue<std::pair<float, labeltype>>> answers(qsize);
    return load_calculated_gt(answers, loc_save_gt, qsize);
}

std::tuple<std::vector<size_t>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> *
test_vs_recall_python(unsigned char *massQ, size_t vecsize, size_t qsize, PQbinDoubleGraph<float> *appr_alg, size_t vecdim,
                      vector<std::priority_queue<std::pair<float, labeltype>>> &answers, size_t k, bool similarity, bool knn_cutoff_applied, float knn_cutoff, std::string index_file, bool hamming, bool weighted, std::string query_type)
{
    return test_vs_recall(massQ, vecsize, qsize, *appr_alg, vecdim, answers, k, similarity, knn_cutoff_applied, knn_cutoff, index_file, hamming, weighted, query_type);
}

static std::tuple<std::vector<size_t>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> *
test_python(PQbinDoubleGraph<float> *appr_alg, py::array_t<float> query_array, bool load_gt, string loc_save_gt, std::string prefix, std::string query_type = "LIGS", int k = 1, string dist_mult = "euclidean")
{

    auto shape = query_array.shape();
    size_t qsize = shape[0];  // Number of rows
    size_t qdim = shape[1];   // Number of rows
    size_t vecdim = shape[1]; // Number of rowspy::array_t<float> query_array = /* Assume this is already defined */;

    // Allocate memory for massQ
    unsigned char* massQ = new unsigned char[qsize * qdim * sizeof(float)];

    // Convert the query array to massQ
    convert_to_queries(query_array, massQ);
    std::vector<std::priority_queue<std::pair<float, labeltype>>> answers(qsize);
    size_t vecsize = appr_alg->max_elements_;
    if (load_gt)
    {
        answers = load_gt_python(loc_save_gt, qsize);
    }
    else
    {
        answers = calculate_gt_python(appr_alg, query_array, loc_save_gt, dist_mult);
    }
    auto results_ptr = test_vs_recall_python(massQ, vecsize, qsize, appr_alg, vecdim, answers, k, false, false, false, prefix, false, false, query_type);
    // Dereference the tuple pointer to access its elements
    auto& results = *results_ptr;

    // Save results to file


    return results_ptr;
}

unsigned char *convert_to_queries_wrapper(py::array_t<float> input_array, unsigned char *massQ)
{
    convert_to_queries(input_array, massQ);
}

PYBIND11_MODULE(python_interface, m)
{
    m.doc() = "Python interface for PQbinDoubleGraph operations, providing functionality for graph initialization, building, node insertion, and deletion, as well as utilities for converting NumPy arrays to C++ data structures.";

    m.def("initialize_graph_python", &initialize_graph_python,
          "Initializes the PQbinDoubleGraph with specified parameters.\n\n"
          "Parameters:\n"
          "    vecsize (int): The size of the vector.\n"
          "    vecdim (int): The dimension of the vector.\n"
          "    subdivision (int): The subdivision parameter for the graph.\n"
          "    centers (int): The number of centers in the graph.\n"
          "    hyperplane_count (int): The number of hyperplances in the graph.\n"
          "    pq_table_no (int): The number of PQ tables.\n"
          "    hyperplane_no (int): The number of hyperplanes.\n"
          "    dist_metric (str): The distance metric to be used ('euclidean', 'angular', 'cosine', 'inner_product').\n\n"
          "Returns:\n"
          "    PQbinDoubleGraph<float>*: A pointer to the initialized graph.",
          py::arg("vecsize"), py::arg("vecdim"), py::arg("subdivision"),
          py::arg("centers"), py::arg("hyperplane_count"), py::arg("pq_table_no"), py::arg("hyperplane_no"),
          py::arg("dist_metric"));

    m.def("build_initial_python", &build_initial_python,
          "Builds the initial graph structure using the provided data.\n\n"
          "Parameters:\n"
          "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
          "    input_array (numpy.ndarray): The input data as a NumPy array.\n"
          "    label (list[int]): The list of labels for the input data.\n"
          "    vecsize (int): The size of the vector.\n"
          "    vecdim (int): The dimension of the vector.\n"
          "    dist_mult (str, optional): The distance metric multiplier ('euclidean' by default).\n"
          "    learn_data (str, optional): Path to learning data (empty string by default).\n"
          "    rand_centroids (bool, optional): Whether to use random centroids (True by default).\n\n"
          "Returns:\n"
          "    float: The build time of the graph.",
          py::arg("appr_alg"), py::arg("input_array"), py::arg("label"),
          py::arg("vecsize"), py::arg("vecdim"), py::arg("dist_mult") = "euclidean",
          py::arg("learn_data") = "", py::arg("rand_centroids") = true);

    m.def("delete_nodes_python", &delete_nodes_python,
          "Deletes nodes from the graph based on the provided list of nodes to delete.\n\n"
          "Parameters:\n"
          "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
          "    nodes_to_delete (list[int]): The list of node indices to delete.\n\n"
          "Returns:\n"
          "    float: The total time taken to delete the nodes.",
          py::arg("appr_alg"), py::arg("nodes_to_delete"));

    m.def("insert_nodes_python", &insert_nodes_python,
          "Inserts new nodes into the graph.\n\n"
          "Parameters:\n"
          "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
          "    input_array (numpy.ndarray): The input data as a NumPy array.\n"
          "    label (list[int]): The list of labels for the input data.\n"
          "    vecdim (int): The dimension of the vector.\n"
          "    dist_mult (str, optional): The distance metric multiplier ('euclidean' by default).\n\n"
          "Returns:\n"
          "    float: The build time of the graph after inserting nodes.",
          py::arg("appr_alg"), py::arg("input_array"), py::arg("label"),
          py::arg("vecdim"), py::arg("dist_mult") = "euclidean");

    m.def("convert_numpy_matrix", &convert_numpy_matrix,
          "Converts a NumPy matrix to a vector of byte arrays, where each byte array represents the bytes of a row in the NumPy matrix.\n\n"
          "Parameters:\n"
          "    input (numpy.ndarray): The input NumPy matrix.\n\n"
          "Returns:\n"
          "    list of byte arrays: A list where each element is a byte array representing a row of the input matrix.",
          py::arg("input"));

    // Continuing from the previous module definition...

    m.def("convert_to_queries", &convert_to_queries_wrapper,
          "Converts a 2D NumPy array of floats to a flat array of queries represented as unsigned char pointers. Each float is reinterpreted as 4 unsigned chars.\n\n"
          "Parameters:\n"
          "    input_array (numpy.ndarray): The input 2D NumPy array of floats.\n"
          "    list of bytes: A list where each element is a bytes object representing a row in the input array, reinterpreted as unsigned chars.",
          py::arg("input_array"), py::arg("massQ"));

    // m.def(
    //     "calculate_gt_python", [](PQbinDoubleGraph<float> *appr_alg, py::array_t<float> query_array, string loc_save_gt, string dist_metric)
    //     { return calculate_gt_python(appr_alg, query_array, loc_save_gt, dist_metric); },
    //     "Calculates and saves the ground truth for a given set of queries using the specified distance metric.\n\n"
    //     "Parameters:\n"
    //     "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
    //     "    query_array (numpy.ndarray): The query data as a NumPy array.\n"
    //     "    loc_save_gt (str): The location (path) where the ground truth should be saved.\n"
    //     "    dist_metric (str): The distance metric to use ('euclidean', 'angular', 'cosine', 'inner_product').\n\n"
    //     "Returns:\n"
    //     "    list of priority queues: A list where each priority queue contains the ground truth (distance, label) pairs for each query.",
    //     py::arg("appr_alg"), py::arg("query_array"), py::arg("loc_save_gt"), py::arg("dist_metric"));

    m.def("load_gt_python", &load_gt_python,
          "Loads the ground truth data from a specified location.\n\n"
          "Parameters:\n"
          "    loc_save_gt (str): The location (path) where the ground truth is saved.\n"
          "    qsize (int): The number of queries (to ensure the loaded ground truth matches the expected query size).\n\n"
          "Returns:\n"
          "    list of priority queues: A list where each priority queue contains the ground truth (distance, label) pairs for each query.",
          py::arg("loc_save_gt"), py::arg("qsize"));

        m.def("check_label", &check_label,
          "checks if label exists.\n\n"
          "Parameters:\n"
          "    appr_alg (str): Graph.\n"
          "    label (int): Label\n\n"
          "Returns:\n"
          "    list of priority queues: A list where each priority queue contains the ground truth (distance, label) pairs for each query.",
          py::arg("appr_alg"), py::arg("label"));

    m.def("get_value", &get_value,
          "fetch the values of the label.\n\n"
          "Parameters:\n"
          "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
          "    nlabel.\n\n"
          "Returns:\n"
          "    float: The total time taken to delete the nodes.",
          py::arg("appr_alg"), py::arg("label"));

    m.def(
        "test_vs_recall_python", [](unsigned char *massQ, size_t vecsize, size_t qsize, PQbinDoubleGraph<float> *appr_alg, size_t vecdim, vector<std::priority_queue<std::pair<float, labeltype>>> &answers, size_t k, bool similarity, bool knn_cutoff_applied, float knn_cutoff, std::string index_file, bool hamming, bool weighted, std::string query_type)
        { return test_vs_recall_python(massQ, vecsize, qsize, appr_alg, vecdim, answers, k, similarity, knn_cutoff_applied, knn_cutoff, index_file, hamming, weighted, query_type); },
        "Tests the recall versus precision for the given set of queries and their ground truths.\n\n"
        "Parameters:\n"
        "    massQ (list of bytes): The query vectors converted to byte format.\n"
        "    vecsize (int): The number of vectors.\n"
        "    qsize (int): The number of queries.\n"
        "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
        "    vecdim (int): The dimension of each vector.\n"
        "    answers (list of priority queues): The ground truth answers for each query.\n"
        "    k (int): The number of nearest neighbors to consider.\n"
        "    similarity (bool): Whether to use similarity instead of distance.\n"
        "    knn_cutoff_applied (bool): Whether a cutoff was applied to k-NN search.\n"
        "    knn_cutoff (float): The cutoff value for k-NN search.\n"
        "    index_file (str): The path to the index file for saving.\n"
        "    hamming (bool): Whether to use Hamming distance for evaluation.\n"
        "    weighted (bool): Whether the graph is weighted.\n"
        "    query_type (str): The type of query.\n\n"
        "Returns:\n"
        "    A tuple containing vectors of test results: sizes, recall, precision, and time.",
        py::arg("massQ"), py::arg("vecsize"), py::arg("qsize"), py::arg("appr_alg"), py::arg("vecdim"), py::arg("answers"),
        py::arg("k"), py::arg("similarity"), py::arg("knn_cutoff_applied"), py::arg("knn_cutoff"), py::arg("index_file"), py::arg("hamming"), py::arg("weighted"), py::arg("query_type"));

    m.def(
        "test_python", [](PQbinDoubleGraph<float> *appr_alg, py::array_t<float> query_array, bool load_gt, string loc_save_gt , std::string prefix, std::string query_type, int k, string dist_metric)
        {     auto results = test_python(appr_alg, query_array, load_gt, loc_save_gt, prefix, query_type, k, dist_metric);
        py::dict py_results;
        py_results["recall"] = py::cast(std::get<1>(*results));
        py_results["precision"] = py::cast(std::get<2>(*results));// Assuming there's an additional metric, change as needed
        delete results; // Make sure to delete the results if they are dynamically allocated
        return py_results; },
        "Performs a test for the given queries using the specified graph and metrics, potentially comparing against loaded ground truth.\n\n"
        "Parameters:\n"
        "    appr_alg (PQbinDoubleGraph<float>*): The graph algorithm object.\n"
        "    query_array (numpy.ndarray): The query data as a NumPy array.\n"
        "    load_gt (bool): Whether to load ground truth from a file.\n"
        "    loc_save_gt (str): The location (path) to save or load the ground truth data.\n"
        "    prefix (str): The prefix for saving the test results.\n"
        "    query_type (str, optional): The type of query ('LIGS' by default).\n"
        "    dist_metric (str, optional): The distance metric ('euclidean' by default).\n"
        "    k (int, optional): The number of nearest neighbors to consider (1 by default).\n\n"
        "Returns:\n"
        "    A tuple containing vectors of test results: sizes, recall, precision, and time.",
        py::arg("appr_alg"), py::arg("query_array"), py::arg("load_gt"), py::arg("loc_save_gt"),
        py::arg("prefix"), py::arg("query_type") = "LIGS", py::arg("k") = 1, py::arg("dist_metric") = "euclidean");


    py::class_<hnswlib::PQbinDoubleGraph<float>, std::shared_ptr<hnswlib::PQbinDoubleGraph<float>>>(m, "PQbinDoubleGraph")
        .def(py::init<hnswlib::SpaceInterface<float>*, hnswlib::SpaceInterface<float>*, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t>(),
             py::arg("s"),
             py::arg("s_sub"),
             py::arg("max_elements"),
             py::arg("vecdim"),
             py::arg("divisions"),
             py::arg("centroids"),
             py::arg("hyperplane_count"),
             py::arg("pq_table_no") = 0,
             py::arg("hyper_table_no") = 0,
             py::arg("M") = 16,
             py::arg("ef") = 32)
        ;
}