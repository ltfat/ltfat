#pragma once

#include <string>
#include <memory>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
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
           unifile.reset( sf_open(file.c_str(), SFM_READ, &fileInfo));
           if (!unifile)
                throw std::runtime_error("Could not open file");
           sf_command (unifile.get(), SFC_SET_NORM_FLOAT, NULL, SF_FALSE);
           sf_command (unifile.get(), SFC_SET_NORM_DOUBLE, NULL, SF_FALSE);
        }

        size_t getReadPos() { return sf_seek  (unifile.get(), 0, SEEK_CUR); }
        size_t getNumSamples() { return fileInfo.frames; }
        int getNumChannels() { return fileInfo.channels; }
        int getSampleRate() { return fileInfo.samplerate; }
        int getFormat(){return fileInfo.format;}

        size_t readSamples(vector<vector<SAMPLE>>& v, size_t numSamplesToRead = 0)
        {
            int reqChannels = v.size();
            size_t reqSamples = max_element(v.begin(),v.end(),
                    [](const auto& v1,const auto& v2){ return v1.size() < v2.size();})[0].size();

            if (numSamplesToRead > 0)
                reqSamples = min(reqSamples,numSamplesToRead);

            size_t samplesRead = 0;
            int chNo = getNumChannels();
            if(buffer.size() != reqSamples*chNo)
                buffer.resize(reqSamples*chNo);

            samplesRead = fwd_readSamples(unifile.get(), buffer.data(), reqSamples*chNo);

            for(int ch = 0; ch < min(chNo,reqChannels); ch++)
                for(size_t l = 0; l < min((size_t)samplesRead,v[ch].size()); l++)
                    v[ch][l] = buffer[chNo*l+ch];

            return samplesRead;
        }

        WavReader(const string& file){ attachFile(file);}
    private:
        unique_ptr<SNDFILE,void(*)(SNDFILE*)> unifile{ nullptr,[](auto* p){ sf_close(p);}};
        SF_INFO fileInfo;
        vector<SAMPLE> buffer;

    // Make non-copyable
    WavReader(const WavReader& a) = delete;
    WavReader& operator=(const WavReader& a) = delete;
    static_assert( is_floating_point<SAMPLE>::value,
                   "Only double and float can be used as sample datatype");
};

sf_count_t fwd_writeSamples(SNDFILE *sndfile, double *ptr, sf_count_t items)
{
   return sf_write_double(sndfile,ptr,items);
}

sf_count_t fwd_writeSamples(SNDFILE *sndfile, float *ptr, sf_count_t items)
{
   return sf_write_float(sndfile,ptr,items);
}

template <typename SAMPLE>
class WavWriter
{
    public:
        void attachFile(const char* file) { attachFile(string(file)); }
        void attachFile(const string& file)
        {

           unifile.reset( sf_open(file.c_str(), SFM_WRITE, &fileInfo));
           if (!unifile)
                throw std::runtime_error("Could not open file");
           sf_command (unifile.get(), SFC_SET_NORM_FLOAT, NULL, SF_FALSE);
           sf_command (unifile.get(), SFC_SET_NORM_DOUBLE, NULL, SF_FALSE);
        }

        size_t getWritePos() { return sf_seek  (unifile.get(), 0, SEEK_CUR); }
        size_t getNumSamples() { return fileInfo.frames; }
        int getNumChannels() { return fileInfo.channels; }
        int getSampleRate() { return fileInfo.samplerate; }

        size_t writeSamples(vector<vector<SAMPLE>>& v, size_t numSamplesToWrite = 0)
        {
            int reqChannels = v.size();
            size_t reqSamples = max_element(v.begin(),v.end(),
                    [](const auto& v1,const auto& v2){ return v1.size() < v2.size();})[0].size();

            if (numSamplesToWrite > 0)
                reqSamples = min(reqSamples,numSamplesToWrite);

            size_t writtenSamples = 0;
            int chNo = getNumChannels();
            if(buffer.size() != reqSamples*chNo)
                buffer.resize(reqSamples*chNo);

            for(int ch = 0; ch < min(chNo,reqChannels); ch++)
                for(size_t l = 0; l < v[ch].size(); l++)
                   buffer[chNo*l+ch] = v[ch][l];

            writtenSamples = fwd_writeSamples(unifile.get(), buffer.data(), reqSamples*chNo);

            return writtenSamples;
        }

        WavWriter(const string& file, int samplerate, int channels,
                  int format=SF_FORMAT_WAV|SF_FORMAT_PCM_16 ):
                  fileInfo{0,samplerate,channels,format,0,0}
                  { attachFile(file);}
    private:
        unique_ptr<SNDFILE,void(*)(SNDFILE*)> unifile{ nullptr,[](auto* p){ sf_close(p);}};
        SF_INFO fileInfo;
        vector<SAMPLE> buffer;

    // Make non-copyable
    WavWriter(const WavWriter& a) = delete;
    WavWriter& operator=(const WavWriter& a) = delete;
    static_assert( is_floating_point<SAMPLE>::value,
                   "Only double and float can be used as sample datatype");
};



