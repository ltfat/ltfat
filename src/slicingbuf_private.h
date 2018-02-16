struct LTFAT_NAME(slicing_processor_state)
{
    LTFAT_NAME(slicing_processor_callback)*
    processorCallback; //!< Custom processor callback
    void* userdata; //!< Callback data
    LTFAT_NAME(fifo_processor_state)* fifo_processor;
    LTFAT_REAL* bufIn;
    LTFAT_REAL* bufOut;
    LTFAT_REAL* bufIn_start;
    LTFAT_REAL* bufOut_start;
    ltfat_int  winLen;
    ltfat_int taperLen;
    ltfat_int zpadLen;
    LTFAT_REAL* ga;
    LTFAT_REAL* gs;

};

int
LTFAT_NAME(slicing_processor_execute_callback)(void* userdata,
        const LTFAT_REAL in[], int winLen, int W, LTFAT_REAL out[]);
