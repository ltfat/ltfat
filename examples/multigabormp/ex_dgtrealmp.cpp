#include "../ltfathelper.h"
#include <sndfile.h>
#include <valarray>
#include <vector>
#include <memory>
using namespace std;

int loadwavfile(const char* name, vector<LTFAT_REAL>& f, int* Ls, int* W)
{
    SF_INFO info;
    SNDFILE* sf = sf_open(name, SFM_READ, &info);
    if (!sf) {
        printf("Failed to load %s\n", name);
        exit(1);
    }

    int maxL = 6*info.samplerate;

    

    *Ls = min(maxL,(int)info.frames) * info.channels;
    *W = info.channels;
    f.resize(*Ls);
#ifdef LTFAT_DOUBLE
    *Ls = sf_read_double(sf, f.data(), *Ls);
#else
    *Ls = sf_read_float(sf, f.data(), *Ls);
#endif
    *Ls /= *W;

    for (int l = 0; l < *Ls; l++)
       f[l] = f[(*W) * l];

    sf_close(sf);
    return 0;
}

char fileTemplate[] = "/home/susnak/Desktop/SQAM/%02d.wav";
/* char file[1024]; */

// const char* file = "/home/susnak/dev/ltfat/signals/gspi.wav";
const char* file = "/home/susnak/Desktop/TSM/Cartoon4s.wav";

int main(int argc, char* argv[])
{
    if (!ltfat_int_is_compatible(sizeof(int)))
    {
        std::cout << "Incompatible size of int. libltfat was probably compiled with -DLTFAT_LARGEARRAYS" << std::endl;
        exit(-1);
    }
    int Ls,W;
    vector<LTFAT_REAL> f;
    vector<unique_ptr<LTFAT_COMPLEX>> coef;

    for (int fileNo = 1; fileNo <= 1; fileNo++) {

        /* snprintf(file, 1023, fileTemplate, fileNo); */
        printf(file);
        printf("\n");
        loadwavfile(file, f, &Ls, &W);
        // printf("Ls=%d,W=%d\n", Ls, W);

        LTFAT_NAME(dgtrealmp_state)*  plan = NULL;
        LTFAT_NAME(dgtrealmp_parbuf)* pbuf = NULL;
        LTFAT_NAME(dgtrealmp_parbuf_init)(&pbuf);

        // LTFAT_NAME(dgtrealmp_setparbuf_alg)(pbuf, ltfat_dgtmp_alg_LocOMP);

        // LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, LTFAT_BLACKMAN, 8192, 2048, 8192);
        // LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, LTFAT_BLACKMAN, 4096, 1024, 4096);
         LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, LTFAT_BLACKMAN, 2048,  512, 2048);
        // LTFAT_NAME(dgtrealmp_parbuf_mod_makelasttight)(pbuf);
        // LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, LTFAT_BLACKMAN, 1024,  256, 1024);
        //  LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, LTFAT_BLACKMAN,  512,  128,  512);
        // LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, LTFAT_BLACKMAN,  256,  64,  256);

        ltfat_int L = LTFAT_NAME(dgtrealmp_parbuf_nextcompatlen)(pbuf,Ls);
        f.resize(L,0.0);
        vector<LTFAT_REAL> fout(L);

        for(int pidx=0;pidx < LTFAT_NAME(dgtrealmp_getparbuf_dictno)(pbuf); pidx++ )
        {
            ltfat_int clen = LTFAT_NAME(dgtrealmp_parbuf_nextcoefsize)(pbuf,Ls,pidx);
            coef.push_back( static_cast<unique_ptr<LTFAT_COMPLEX>>(LTFAT_NAME_COMPLEX(malloc)(clen)) );
        }

        LTFAT_NAME(dgtrealmp_setparbuf_phaseconv)(pbuf, LTFAT_TIMEINV);
        LTFAT_NAME(dgtrealmp_setparbuf_snrdb)(pbuf, 50);
        LTFAT_NAME(dgtrealmp_setparbuf_kernrelthr)(pbuf, 1e-4);
        LTFAT_NAME(dgtrealmp_setparbuf_maxatoms)(pbuf, 0.8*L);
        LTFAT_NAME(dgtrealmp_setparbuf_iterstep)(pbuf, L);

        LTFAT_NAME(dgtrealmp_init)( pbuf, L, &plan);

        LTFAT_NAME(dgtrealmp_parbuf_done)(&pbuf);

        auto t1 = Clock::now();
        LTFAT_NAME(dgtrealmp_execute)(plan, f.data(), (LTFAT_COMPLEX**) coef.data(), fout.data());
        auto t2 = Clock::now();
        int dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        cout << "DURATION: " << dur << " ms" << std::endl;

        size_t atoms;
        LTFAT_NAME(dgtrealmp_get_numatoms)(plan,&atoms);
        size_t iters;
        LTFAT_NAME(dgtrealmp_get_numiters)(plan,&iters);

        LTFAT_NAME(dgtrealmp_done)(&plan);

        LTFAT_REAL snr;
        LTFAT_NAME(snr)(f.data(), fout.data(), L, &snr);

        printf("atoms=%u,iters=%u,SNR=%2.3f dB, perit=%2.3f us\n",atoms, iters, snr, 1000*dur/((double)iters) );
    }

    return 0;
}
