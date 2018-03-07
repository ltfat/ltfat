#pragma once

#include <string>
#include <memory>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <sndfile.h>
using namespace std;

sf_count_t fwd_readSamples(SNDFILE *sndfile, double *ptr, sf_count_t items)
{
   return sf_read_double(sndfile,ptr,items);
}

sf_count_t fwd_readSamples(SNDFILE *sndfile, float *ptr, sf_count_t items)
{
   return sf_read_float(sndfile,ptr,items);
}

template <typename SAMPLE>
class WavReader
{
    public:
        void attachFile(const char* file) { attachFile(string(file)); }
        void attachFile(const string& file)
        {
            unifile = unique_ptr<SNDFILE,void(*)(SNDFILE*)>(
                    sf_open(file.c_str(), SFM_READ, &fileInfo),
                    [](SNDFILE* p){ sf_close(p);});
           if (!unifile)
                throw std::runtime_error("Could not open file");
        }

        size_t getReadPos() { return sf_seek  (unifile.get(), 0, SEEK_CUR); }
        size_t getNumSamples() { return fileInfo.frames; }
        int getNumChannels() { return fileInfo.channels; }
        int getSampleRate() { return fileInfo.samplerate; }

        size_t readSamples(vector<SAMPLE>& v)
        {
            size_t toret = 0;
            int chNo = getNumChannels();
            if(buffer.size() != v.size()*chNo)
                buffer.resize(v.size()*chNo);

            toret = fwd_readSamples(unifile.get(), buffer.data(), v.size()*chNo);

            for(int ch = 0; ch < chNo; ch++)
            {
                for(size_t l = 0; l < toret; l++)
                {
                    v[l] = buffer[chNo*l];
                }
            }

            return toret;
        }

        WavReader(const string& file){ attachFile(file);}
    private:
        unique_ptr<SNDFILE,void(*)(SNDFILE*)> unifile{ nullptr, nullptr };
        SF_INFO fileInfo;
        vector<SAMPLE> buffer;

    // Make non-copyable
    WavReader(const WavReader& a) = delete;
    WavReader& operator=(const WavReader& a) = delete;
    // static_assert( typeid(SAMPLE) == typeid(double) ||
    //                typeid(SAMPLE) == typeid(float),
    //                "Only double and float can be used as sample datatype");
};



