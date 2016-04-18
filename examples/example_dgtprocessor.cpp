#include "ltfat.h"

class DGTProcessor {
public:
    class Callback {
    public:
        virtual int callback(const std::complex<float> in[], const int M2,
            std::complex<float> out[]) noexcept = 0;
    };

    DGTProcessor(LTFAT_FIRWIN win, int gl, int a, int M)
        : g(win)
        , gl(gl)
        , a(a)
        , M(M)
    {
    }

    virtual ~DGTProcessor()
    {
        if (processor_struct)
            rtdgtreal_processor_done_s(processor_struct);
        if (callback)
            delete callback;
    }

    void registerCallback(DGTProcessor::Callback* callback)
    {
        this->callback = callback;
        if (callback) {
            if (processor_struct)
                rtdgtreal_processor_done_s(processor_struct);
            processor_struct = rtdgtreal_processor_wininit_s(
                g, gl, a, M, DGTProcessor::callbackWrapper, static_cast<void*>(callback));
        }
    }

    int process(const float in[], int bufLen, float out[]) noexcept
    {
        if (processor_struct)
            rtdgtreal_processor_execute_s(processor_struct, in, bufLen, out);
        else
            memcpy(out, in, bufLen * sizeof(float));
    }

private:
    rtdgtreal_processor_s* processor_struct{ nullptr };
    DGTProcessor::Callback* callback{ nullptr };
    LTFAT_FIRWIN g;
    int gl;
    int a;
    int M;

    static void callbackWrapper(
        void* userdata, const float in[][2], const int M2, float out[][2])
    {
        static_cast<DGTProcessor::Callback*>(userdata)->callback(
            reinterpret_cast<const std::complex<float>*>(in), M2,
            reinterpret_cast<std::complex<float>*>(out));
    }
};

class SimpleDGTProcessorCallback : public DGTProcessor::Callback {
public:
    int callback(
        const std::complex<float> in[], const int M2, std::complex<float> out) noexcept
    {
    }
};
