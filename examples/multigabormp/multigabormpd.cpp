#include "ltfathelper.h"
#include "cxxopts.hpp"
#include "wavhandler.h"


template<class T>
using uni_ptrdel = unique_ptr<T, void(*)( T*)>;

int main(int argc, char* argv[])
{
    if (!ltfat_int_is_compatible(sizeof(int)))
    {
        std::cout << "Incompatible size of int. libltfat was probably"
                     " compiled with -DLTFAT_LARGEARRAYS" << std::endl;
        exit(1);
    }

    string inFile, outFile, resFile;
    double targetsnrdb = 40;
    size_t maxit = 0, maxat = 0;
    double seglen = 10;
    double kernthr = 1e-4;
    vector<tuple<string,LTFAT_FIRWIN,int,int>> dicts;
    size_t numSamples = 0;
    int numChannels = 0;
    int sampRate = 0;
    bool do_pedanticsearch = false;

    try
    {
        string examplestr{"Usage:\n" + 
            string(argv[0]) + " -s snr -d win,a,M input.wav"
            + "\nExample:\n" +
            string(argv[0]) + " -s 40 -d HANN,512,2048 input.wav"
        };
        cxxopts::Options options(argv[0], " - example command line options");
        options
        .positional_help("-s snr -d win,a,M input.wav")
        .show_positional_help();

        options.add_options()
        ("i,input", "Input *.wav file", cxxopts::value<string>())
        ("o,output","Output *.wav file", cxxopts::value<string>())
        ("r,residual","Residual *.wav file", cxxopts::value<string>() )
        ("d,dict","Dictionary specification. Format: win1,hop1,channels1:win2,hop2,channels2. "
                  "Example: blackman,512,2048:blackman,256,1024:... Supported windows are: blackman, hann",
                  cxxopts::value<string>() )
        ("s,snr", "Target signal-to-noise ratio",
         cxxopts::value<double>()->default_value(to_string(targetsnrdb)))
        ("maxit", "Maximum number of iterations", cxxopts::value<size_t>() )
        ("maxat", "Maximum number of atoms", cxxopts::value<size_t>() )
        ("kernthr", "Kernel truncation threshold",
         cxxopts::value<double>()->default_value(to_string(kernthr)))
        ("seglen", "Segment length in seconds. 0 disables the segmentation.",
         cxxopts::value<double>()->default_value(to_string(seglen)) )
        ("pedanticsearch", "Enables pedantic search.",
         cxxopts::value<bool>(do_pedanticsearch) )
        ("help", "Print help");

        options.parse_positional({"input"});

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            cout << options.help({""}) << endl;
            exit(0);
        }

        if (result.count("input"))
        {
            inFile = result["input"].as<string>();
            try
            {
                WavReader<LTFAT_REAL> wrtmp{inFile};
                numSamples = wrtmp.getNumSamples();
                numChannels = wrtmp.getNumChannels();
                sampRate = wrtmp.getSampleRate();
            }
            catch(...)
            {
                cout << "Cannot open " << inFile << endl;
                exit(1);
            }
        }
        else
        {
            cout << "No input file specified." << endl;
            cout << examplestr << endl;
            exit(1);
        }

        if (result.count("output"))
            outFile = result["output"].as<string>();

        if (result.count("residual"))
            resFile = result["residual"].as<string>();

        if(result.count("seglen"))
        {
            seglen =  result["seglen"].as<double>();
        }

        if (result.count("snr"))
            targetsnrdb = result["snr"].as<double>();
        else
        {
            cout << "No taget SNR specified." << endl;
            cout << examplestr << endl;
            exit(1);
        }

        if (result.count("dict"))
        {
             string toparse = result["dict"].as<string>() + ":";

             int pos;
             while ((pos = toparse.find(":")) != -1)
             {
                 string dictstr = toparse.substr(0,pos);
                 toparse = toparse.substr(pos+1,toparse.size()-pos);
                 if( !dictstr.empty() )
                 {
                    dictstr += ",";
                    vector<string> dictvec;
                    int pos2;
                    while ((pos2 = dictstr.find(",")) != -1)
                    {
                        string itemstr = dictstr.substr(0,pos2);
                        dictstr = dictstr.substr(pos2+1,dictstr.size()-pos2);
                        if( !itemstr.empty())
                        {
                            dictvec.push_back(itemstr);
                            cout << itemstr << endl;
                        }
                    }
                    if(dictvec.size() != 3)
                    {
                        cout << "Parse error: Dictionary should consist of 3 items: win,a,M" << endl;
                        exit(1);
                    }
                    dicts.push_back(make_tuple(dictvec[0],LTFAT_BLACKMAN,stoi(dictvec[1]),stoi(dictvec[2])));
                 }
             }
        }
        else
        {
            cout << "No dictionary specified." << endl;
            cout << examplestr << endl;
            exit(1);
        }

        if (!result.count("maxat"))
        {
            maxat = (size_t) (numSamples);
        }

        if (!result.count("maxit"))
        {
            maxit = (size_t) (0.8*numSamples);
            maxat = maxit;
        }

        if (result.count("kernthr"))
        {
            kernthr = result["kernthr"].as<double>();
        }
    }
    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }


    LTFAT_NAME(dgtrealmp_parbuf)* pbuf = NULL;
    LTFAT_NAME(dgtrealmp_parbuf_init)(&pbuf);
    auto unipb = uni_ptrdel<LTFAT_NAME(dgtrealmp_parbuf)>(
    pbuf,[](auto* p){ LTFAT_NAME(dgtrealmp_parbuf_done)(&p); });

    // LTFAT_NAME(dgtrealmp_setparbuf_alg)(pbuf, ltfat_dgtmp_alg_LocOMP);

    for(auto dict:dicts)
    {
       if( 0 > LTFAT_NAME(dgtrealmp_parbuf_add_firwin)(pbuf, get<1>(dict), get<3>(dict), get<2>(dict), get<3>(dict)))
       {
          cout << "Bad dictionary: " << get<0>(dict) << "," << get<2>(dict)  << "," << get<3>(dict) << endl;
          exit(1);
       }
    }

    cout << "seglen: " << seglen << endl;
    if( numSamples <= seglen*sampRate )
    {
    ltfat_int L = LTFAT_NAME(dgtrealmp_getparbuf_siglen)(pbuf, numSamples);

    vector<LTFAT_REAL> f(L,0.0);
    vector<LTFAT_REAL> fout(L);
    vector<unique_ptr<LTFAT_COMPLEX[]>> coef;

    for (int pidx = 0; pidx < LTFAT_NAME(dgtrealmp_getparbuf_dictno)(pbuf); pidx++ )
    {
        ltfat_int clen = LTFAT_NAME(dgtrealmp_getparbuf_coeflen)(pbuf, numSamples, pidx);
        coef.push_back( unique_ptr<LTFAT_COMPLEX[]>(new LTFAT_COMPLEX[clen]) );
    }

    WavReader<LTFAT_REAL> wr{inFile};
    wr.readSamples(f);

    LTFAT_NAME(dgtrealmp_setparbuf_phaseconv)(pbuf, LTFAT_TIMEINV);
    LTFAT_NAME(dgtrealmp_setparbuf_pedanticsearch)(pbuf, do_pedanticsearch);
    LTFAT_NAME(dgtrealmp_setparbuf_snrdb)(pbuf, targetsnrdb);
    LTFAT_NAME(dgtrealmp_setparbuf_kernrelthr)(pbuf, kernthr);
    LTFAT_NAME(dgtrealmp_setparbuf_maxatoms)(pbuf, 0.8 * L);
    LTFAT_NAME(dgtrealmp_setparbuf_maxit)(pbuf, 2e5);
    LTFAT_NAME(dgtrealmp_setparbuf_iterstep)(pbuf, L);

    LTFAT_NAME(dgtrealmp_state)*  plan = NULL;
    auto t1 = Clock::now();
    LTFAT_NAME(dgtrealmp_init)( pbuf, L, &plan);
    auto t2 = Clock::now();
    int dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    cout << "INIT DURATION: " << dur << " ms" << std::endl;

    t1 = Clock::now();
    LTFAT_NAME(dgtrealmp_execute_decompose)(plan, f.data(), (LTFAT_COMPLEX**) coef.data());
    t2 = Clock::now();
    dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    cout << "DURATION: " << dur << " ms" << std::endl;

    t1 = Clock::now();
    LTFAT_NAME(dgtrealmp_execute_synthesize)(plan, (const LTFAT_COMPLEX**) coef.data(), NULL, fout.data());
    t2 = Clock::now();
    int dur2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    cout << "SYN DURATION: " << dur2 << " ms" << std::endl;


    size_t atoms;
    LTFAT_NAME(dgtrealmp_get_numatoms)(plan, &atoms);
    size_t iters;
    LTFAT_NAME(dgtrealmp_get_numiters)(plan, &iters);

    LTFAT_NAME(dgtrealmp_done)(&plan);

    LTFAT_REAL snr;
    LTFAT_NAME(snr)(f.data(), fout.data(), L, &snr);

    printf("atoms=%u,iters=%u,SNR=%2.3f dB, perit=%2.3f us\n", atoms, iters, snr, 1000 * dur / ((double)iters) );
    }
    else
    {
    cout << "Unsupported" << endl;
    }
    return 0;
}
