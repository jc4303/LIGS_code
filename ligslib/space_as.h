#pragma once
#include "hnswlib.h"
#include <math.h>  
namespace hnswlib {

    static float
    AngularSim(const void *pVect1, const void *pVect2, const void *qty_ptr) {
        size_t qty = *((size_t *) qty_ptr);
        float res = 0;
        float pTot1 = 0;
        float pTot2 = 0;
        for (unsigned i = 0; i < qty; i++) {
            pTot1 += ((float *) pVect1)[i] *(( float *) pVect1)[i];
            pTot2 += ((float *) pVect2)[i] *(( float *) pVect2)[i];

        }
        pTot1 = sqrt(pTot1);
        pTot2 = sqrt(pTot2);


        for (unsigned i = 0; i < qty; i++) {
            res += ((float *) pVect1)[i] * ((float *) pVect2)[i];
        }
        res = res/pTot2/pTot1;
        return res;

    }

    static float
    AngularDistance(const void *pVect1, const void *pVect2, const void *qty_ptr) {
        return 1.0f - AngularSim(pVect1, pVect2, qty_ptr);
    }

#if defined(USE_AVX)

// Favor using AVX if available.
    static float
    AngularSIMD4ExtAVX(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        float PORTABLE_ALIGN32 TmpRes[8];
        float *pVect1 = (float *) pVect1v;
        float *pVect2 = (float *) pVect2v;
        size_t qty = *((size_t *) qty_ptr);

        size_t qty16 = qty / 16;
        size_t qty4 = qty / 4;

        const float *pEnd1 = pVect1 + 16 * qty16;
        const float *pEnd2 = pVect1 + 4 * qty4;

        __m256 sum256 = _mm256_set1_ps(0);

        while (pVect1 < pEnd1) {
            //_mm_prefetch((char*)(pVect2 + 16), _MM_HINT_T0);

            __m256 v1 = _mm256_loadu_ps(pVect1);
            pVect1 += 8;
            __m256 v2 = _mm256_loadu_ps(pVect2);
            pVect2 += 8;
            sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v1, v2));

            v1 = _mm256_loadu_ps(pVect1);
            pVect1 += 8;
            v2 = _mm256_loadu_ps(pVect2);
            pVect2 += 8;
            sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v1, v2));
        }

        __m128 v1, v2;
        __m128 sum_prod = _mm_add_ps(_mm256_extractf128_ps(sum256, 0), _mm256_extractf128_ps(sum256, 1));

        while (pVect1 < pEnd2) {
            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));
        }

        _mm_store_ps(TmpRes, sum_prod);
        float sum = TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3];;
        return sum;
    }
    
    static float
    AngularDistanceSIMD4ExtAVX(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        return 1.0f - AngularSIMD4ExtAVX(pVect1v, pVect2v, qty_ptr);
    }

#endif

#if defined(USE_SSE)

    static float
    AngularSIMD4ExtSSE(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        float PORTABLE_ALIGN32 TmpRes[8];
        float *pVect1 = (float *) pVect1v;
        float *pVect2 = (float *) pVect2v;
        size_t qty = *((size_t *) qty_ptr);

        size_t qty16 = qty / 16;
        size_t qty4 = qty / 4;

        const float *pEnd1 = pVect1 + 16 * qty16;
        const float *pEnd2 = pVect1 + 4 * qty4;

        __m128 v1, v2;
        __m128 sum_prod = _mm_set1_ps(0);

        while (pVect1 < pEnd1) {
            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));

            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));

            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));

            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));
        }

        while (pVect1 < pEnd2) {
            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));
        }

        _mm_store_ps(TmpRes, sum_prod);
        float sum = TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3];

        return sum;
    }

    static float
    AngularDistanceSIMD4ExtSSE(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        return 1.0f - AngularSim(pVect1v, pVect2v, qty_ptr);
    }

#endif


#if defined(USE_AVX512)

    static float
    AngularSIMD16ExtAVX512(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        float PORTABLE_ALIGN64 TmpRes[16];
        float *pVect1 = (float *) pVect1v;
        float *pVect2 = (float *) pVect2v;
        size_t qty = *((size_t *) qty_ptr);

        size_t qty16 = qty / 16;


        const float *pEnd1 = pVect1 + 16 * qty16;

        __m512 sum512 = _mm512_set1_ps(0);

        while (pVect1 < pEnd1) {
            //_mm_prefetch((char*)(pVect2 + 16), _MM_HINT_T0);

            __m512 v1 = _mm512_loadu_ps(pVect1);
            pVect1 += 16;
            __m512 v2 = _mm512_loadu_ps(pVect2);
            pVect2 += 16;
            sum512 = _mm512_add_ps(sum512, _mm512_mul_ps(v1, v2));
        }

        _mm512_store_ps(TmpRes, sum512);
        float sum = TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3] + TmpRes[4] + TmpRes[5] + TmpRes[6] + TmpRes[7] + TmpRes[8] + TmpRes[9] + TmpRes[10] + TmpRes[11] + TmpRes[12] + TmpRes[13] + TmpRes[14] + TmpRes[15];

        return sum;
    }

    static float
    AngularDistanceSIMD16ExtAVX512(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        return 1.0f - AngularSIMD16ExtAVX512(pVect1v, pVect2v, qty_ptr);
    }

#endif

#if defined(USE_AVX)

    static float
    AngularSIMD16ExtAVX(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        float PORTABLE_ALIGN32 TmpRes[8];
        float *pVect1 = (float *) pVect1v;
        float *pVect2 = (float *) pVect2v;
        size_t qty = *((size_t *) qty_ptr);

        size_t qty16 = qty / 16;


        const float *pEnd1 = pVect1 + 16 * qty16;

        __m256 sum256 = _mm256_set1_ps(0);

        while (pVect1 < pEnd1) {
            //_mm_prefetch((char*)(pVect2 + 16), _MM_HINT_T0);

            __m256 v1 = _mm256_loadu_ps(pVect1);
            pVect1 += 8;
            __m256 v2 = _mm256_loadu_ps(pVect2);
            pVect2 += 8;
            sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v1, v2));

            v1 = _mm256_loadu_ps(pVect1);
            pVect1 += 8;
            v2 = _mm256_loadu_ps(pVect2);
            pVect2 += 8;
            sum256 = _mm256_add_ps(sum256, _mm256_mul_ps(v1, v2));
        }

        _mm256_store_ps(TmpRes, sum256);
        float sum = TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3] + TmpRes[4] + TmpRes[5] + TmpRes[6] + TmpRes[7];

        return sum;
    }

    static float
    AngularDistanceSIMD16ExtAVX(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        return 1.0f - AngularSIMD16ExtAVX(pVect1v, pVect2v, qty_ptr);
    }

#endif

#if defined(USE_SSE)

    static float
    AngularSIMD16ExtSSE(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        float PORTABLE_ALIGN32 TmpRes[8];
        float *pVect1 = (float *) pVect1v;
        float *pVect2 = (float *) pVect2v;
        size_t qty = *((size_t *) qty_ptr);

        size_t qty16 = qty / 16;

        const float *pEnd1 = pVect1 + 16 * qty16;

        __m128 v1, v2;
        __m128 sum_prod = _mm_set1_ps(0);

        while (pVect1 < pEnd1) {
            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));

            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));

            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));

            v1 = _mm_loadu_ps(pVect1);
            pVect1 += 4;
            v2 = _mm_loadu_ps(pVect2);
            pVect2 += 4;
            sum_prod = _mm_add_ps(sum_prod, _mm_mul_ps(v1, v2));
        }
        _mm_store_ps(TmpRes, sum_prod);
        float sum = TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3];

        return sum;
    }

    static float
    AngularDistanceSIMD16ExtSSE(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        return 1.0f - AngularSim(pVect1v, pVect2v, qty_ptr);
    }

#endif

#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_AVX512)
    DISTFUNC<float> AngularSIMD16Ext = AngularSim;
    DISTFUNC<float> AngularSIMD4Ext = AngularSim;
    DISTFUNC<float> AngularDistanceSIMD16Ext = AngularSim;
    DISTFUNC<float> AngularDistanceSIMD4Ext = AngularSim;

    static float
    AngularDistanceSIMD16ExtResiduals(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        size_t qty = *((size_t *) qty_ptr);
        size_t qty16 = qty >> 4 << 4;
        float res = AngularSIMD16Ext(pVect1v, pVect2v, &qty16);
        float *pVect1 = (float *) pVect1v + qty16;
        float *pVect2 = (float *) pVect2v + qty16;

        size_t qty_left = qty - qty16;
        float res_tail = AngularSim(pVect1, pVect2, &qty_left);
        return 1.0f - (res + res_tail);
    }

    static float
   AngularDistanceSIMD4ExtResiduals(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        size_t qty = *((size_t *) qty_ptr);
        size_t qty4 = qty >> 2 << 2;

        float res = AngularSIMD4Ext(pVect1v, pVect2v, &qty4);
        size_t qty_left = qty - qty4;

        float *pVect1 = (float *) pVect1v + qty4;
        float *pVect2 = (float *) pVect2v + qty4;
        float res_tail = AngularSim(pVect1, pVect2, &qty_left);

        return 1.0f - (res + res_tail);
    }
#endif

    class AngularSpace : public SpaceInterface<float> {

        DISTFUNC<float> fstdistfunc_;
        size_t data_size_;
        size_t dim_;
    public:
        AngularSpace(size_t dim) {
            fstdistfunc_ = AngularDistance;
    #if defined(USE_AVX) || defined(USE_SSE) || defined(USE_AVX512)
        #if defined(USE_AVX512)
            if (AVX512Capable()) {
                AngularSIMD16Ext = AngularSIMD16ExtAVX512;
                AngularDistanceSIMD16Ext = AngularDistanceSIMD16ExtAVX512;
            } else if (AVXCapable()) {
                AngularSIMD16Ext = AngularSIMD16ExtAVX;
                AngularDistanceSIMD16Ext = AngularDistanceSIMD16ExtAVX;
            }
        #elif defined(USE_AVX)
            if (AVXCapable()) {
                AngularSIMD16Ext = AngularSIMD16ExtAVX;
                AngularDistanceSIMD16Ext = AngularDistanceSIMD16ExtAVX;
            }
        #endif
        #if defined(USE_AVX)
            if (AVXCapable()) {
                AngularSIMD4Ext = AngularSIMD4ExtAVX;
               AngularDistanceSIMD4Ext = AngularDistanceSIMD4ExtAVX;
            }
        #endif

            if (dim % 16 == 0)
                fstdistfunc_ = AngularDistance;
            else if (dim % 4 == 0)
                fstdistfunc_ = AngularDistance;
            else if (dim > 16)
                fstdistfunc_ = AngularDistance;
            else if (dim > 4)
                fstdistfunc_ = AngularDistance;
    #endif
            dim_ = dim;
            data_size_ = dim * sizeof(float);
        }

        size_t get_data_size() {
            return data_size_;
        }

        DISTFUNC<float> get_dist_func() {
            return fstdistfunc_;
        }

        void *get_dist_func_param() {
            return &dim_;
        }

    ~AngularSpace() {}
    };

}
