#include "../ltfathelper.h"
#include <sndfile.h>

int loadwavfile(const char* name, LTFAT_REAL** f, int* Ls, int* W)
{
    SF_INFO info;
    SNDFILE* sf = sf_open(name, SFM_READ, &info);
    if (!sf) {
        printf("Failed to load %s\n", name);
        exit(1);
    }

    *Ls = info.frames * info.channels;
    *W = info.channels;
    LTFAT_REAL* ftmp = LTFAT_NAME(calloc)(*Ls);
    *f = LTFAT_NAME(calloc)(*Ls / (*W));

#ifdef LTFAT_DOUBLE
    *Ls = sf_read_double(sf, ftmp, *Ls);
#else
    *Ls = sf_read_float(sf, ftmp, *Ls);
#endif
    *Ls /= *W;

        for (int l = 0; l < *Ls; l++) {
            (*f)[l] = ftmp[(*W) * l];
        }

    ltfat_safefree(ftmp);
    sf_close(sf);
    return 0;
}

char fileTemplate[] = "/home/susnak/Desktop/SQAM/%02d.wav";
/* char file[1024]; */
const char* file = "/home/susnak/dev/ltfat/signals/gspi.wav";

int main(int argc, char* argv[])
{
    if (!ltfat_int_is_compatible(sizeof(int)))
    {
        std::cout << "Incompatible size of int. Compile libltfat with -DLTFAT_COMPAT32" << std::endl;
        exit(-1);
    }
    int Ls, L;
    int W;
    LTFAT_REAL* f;
    LTFAT_COMPLEX* cout;
    int fs = 44100;

    for (int fileNo = 1; fileNo <= 1; fileNo++) {

        /* snprintf(file, 1023, fileTemplate, fileNo); */
        printf(file);
        printf("\n");
        loadwavfile(file, &f, &Ls, &W);
        /* printf("Ls=%d,W=%d\n", Ls, W); */

        ltfat_int a = 512;
        ltfat_int M = 2048;
        ltfat_int gl = 2048;
        ltfat_int M2 = M / 2 + 1;

        L = ltfat_dgtlength(Ls, a, M);
        printf("L=%td\n", L);

        f = LTFAT_NAME(postpad)(f, Ls, L);

        LTFAT_REAL* fout = LTFAT_NAME(malloc)(L);

        LTFAT_REAL* g = LTFAT_NAME(malloc)(gl);
        cout = LTFAT_NAME_COMPLEX(malloc)(M2 * L / a);

        LTFAT_NAME(firwin)(LTFAT_BLACKMAN, gl, g);
        LTFAT_NAME(normalize)(g, gl, LTFAT_NORM_ENERGY, g);

        LTFAT_NAME(dgtrealmp_state)* plan = NULL;

        ltfat_dgtmp_params* params = ltfat_dgtmp_params_allocdef();
        ltfat_dgtmp_setpar_phaseconv(params, LTFAT_TIMEINV);
        ltfat_dgtmp_setpar_errtoldb(params, -40);
        ltfat_dgtmp_setpar_kernrelthr(params, 1e-4);
        ltfat_dgtmp_setpar_maxatoms(params, 0.8 * L);
        ltfat_dgtmp_setpar_iterstep(params, 1e6);
        // ltfat_dgtrealmp_setpar_alg(params, ltfat_dgtrealmp_alg_LocCyclicMP);

        int retval = LTFAT_NAME(dgtrealmp_init_gen)(
            (const LTFAT_REAL**)&g, &gl, L, 1, &a, &M, params, &plan);
        printf("status=%d\n", retval);

        ltfat_dgtmp_params_free(params);

        auto t1 = Clock::now();

        retval = LTFAT_NAME(dgtrealmp_execute)(plan, f, &cout, fout);
        auto t2 = Clock::now();
        std::cout
            << "Delta t2-t1: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
            << " miliseconds" << std::endl;
        printf("status=%d\n", retval);

        size_t atoms;
        LTFAT_NAME(dgtrealmp_get_numatoms)(plan,&atoms);
        size_t iters;
        LTFAT_NAME(dgtrealmp_get_numiters)(plan,&iters);

        LTFAT_NAME(dgtrealmp_done)(&plan);

        LTFAT_REAL snr;
        LTFAT_NAME(snr)(f, fout, L, &snr);

        printf("atoms=%u,iters=%u,SNR=%2.3f dB\n",atoms, iters, snr);
        ltfat_safefree(g);
        ltfat_safefree(f);
        ltfat_safefree(cout);
        ltfat_safefree(fout);
    }

    return 0;
}
