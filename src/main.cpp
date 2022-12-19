#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <ipp.h>

const double Umax = 32767.0; //32767;

int main(int argc, char const *argv[])
{
    // std::ofstream lfm_float_re("lfm_re_float.dat", std::ios_base::out | std::ios_base::binary);
    // std::ofstream lfm_float_im("lfm_im_float.dat", std::ios_base::out | std::ios_base::binary);

    // ��� ���������� � Vivado 
    std::ofstream lfm_int16_re("lfm_re_int16.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_int16_im("lfm_im_int16.dat", std::ios_base::out | std::ios_base::binary);

    std::ofstream lfm_int16_shift_re("lfm_int16_shift_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_int16_shift_im("lfm_int16_shift_im.dat", std::ios_base::out | std::ios_base::binary);

    // ��� ������������� ���������� - ���
    std::ofstream lfm_re("lfm_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_im("lfm_im.dat", std::ios_base::out | std::ios_base::binary);

    // ��� ������������� ���������� - ��������� ��� ������
    std::ofstream lfm_shift_re("lfm_shift_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_shift_im("lfm_shift_im.dat", std::ios_base::out | std::ios_base::binary);

    // ��� ������������� ���������� - ��������� ����������
    std::ofstream correl_re("correl_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream correl_im("correl_im.dat", std::ios_base::out | std::ios_base::binary);

    double freqDev = 10000000;
    double samplingFrequency = 100000000;
    double signalDuration = 4 / 1000000.0;
    double fi = 0;

    int size = 10;
    int shift = 256;

    Ipp32fc* data = nullptr;
    Ipp32fc* data_shift = nullptr;
    Ipp32fc* m_bufferFFTRange = nullptr;

    int sizeDFTSpec = 0;
    int sizeDFTInitBuf = 0;
    int sizeDFTWorkBuf = 0;

    int fft_len = std::pow(2, size);

    Ipp8u* pDFTSpec = nullptr;
    Ipp8u* pDFTInitBuf = nullptr;
    Ipp8u* pDFTWorkBuf = nullptr;
    IppsFFTSpec_C_32fc* pSpec = 0;

    ippsFFTGetSize_C_32fc(size, IPP_NODIV_BY_ANY, ippAlgHintAccurate, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);

    pDFTSpec = ippsMalloc_8u(sizeDFTSpec);
    pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
    pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);
    
    ippsFFTInit_C_32fc(&pSpec, size, IPP_NODIV_BY_ANY, ippAlgHintAccurate, pDFTSpec, pDFTInitBuf);

    data = ippsMalloc_32fc(fft_len);
    data_shift = ippsMalloc_32fc(fft_len);

    ippsZero_32fc(data, fft_len);
    ippsZero_32fc(data_shift, fft_len);

    double lengthEtalonSignal = signalDuration * samplingFrequency;     // ����������� ������ ���������� �������
    for (int i = 0; i < lengthEtalonSignal; ++i)                        // ������������ �������
    {
        double samplingT = (double)i / (double)samplingFrequency - (double)signalDuration / 2.0;    // ���������� ������� �������������

        // ������������ �������������� ������������
        double tempRe = round(Umax * cos(2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration));                       //TODO - �������� �� cos 
        if (tempRe > Umax)
            tempRe = Umax;
        else if (tempRe < -Umax)
            tempRe = -Umax;

        // ������������ ������ ������������
        double tempIm = round(Umax * sin(2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration + fi * M_PI / 180.0));   //TODO - �������� �� sin
        if (tempIm > Umax)
            tempIm = Umax;
        else if (tempIm < -Umax)
            tempIm = -Umax;

        data[i] = { (float)tempRe, (float)tempIm };
        data_shift[i + shift] = { (float)tempRe, (float)tempIm };
    }

    for (size_t i = 0; i < fft_len; i++)
    {
        lfm_re << i << '\t' << data[i].re << std::endl;
        lfm_im << i << '\t' << data[i].im << std::endl;
        lfm_shift_re << i << '\t' << data_shift[i].re << std::endl;
        lfm_shift_im << i << '\t' << data_shift[i].im << std::endl;
    }

    for (size_t i = 0; i < fft_len; i++)
    {
        int16_t tempRe_int16 = data[i].re;
        int16_t tempIm_int16 = data[i].im;

        const unsigned char* pf_re = reinterpret_cast<const unsigned char*>(&tempRe_int16);
        const unsigned char* pf_im = reinterpret_cast<const unsigned char*>(&tempIm_int16);

        lfm_int16_re << pf_re[1];
        lfm_int16_re << pf_re[0];

        lfm_int16_im << pf_im[1];
        lfm_int16_im << pf_im[0];
    }

    for (size_t i = 0; i < fft_len; i++)
    {
        int16_t tempRe_int16 = data_shift[i].re;
        int16_t tempIm_int16 = data_shift[i].im;

        const unsigned char* pf_re = reinterpret_cast<const unsigned char*>(&tempRe_int16);
        const unsigned char* pf_im = reinterpret_cast<const unsigned char*>(&tempIm_int16);

        lfm_int16_shift_re << pf_re[1];
        lfm_int16_shift_re << pf_re[0];

        lfm_int16_shift_im << pf_im[1];
        lfm_int16_shift_im << pf_im[0];
    }

    ippsFFTFwd_CToC_32fc(data, data, pSpec, pDFTWorkBuf);

    ippsFFTFwd_CToC_32fc(data_shift, data_shift, pSpec, pDFTWorkBuf);

    Ipp32fc* correl = nullptr;

    correl = ippsMalloc_32fc(fft_len);

    ippsZero_32fc(correl, fft_len);

    ippsConj_32fc(data, data, fft_len);

    /*
    for (size_t i = 0; i < fft_len; i++)
    {
        const unsigned char* float_re = reinterpret_cast<const unsigned char*>(&correl[i].re);
        const unsigned char* float_im = reinterpret_cast<const unsigned char*>(&correl[i].im);

        lfm_float_re << float_re[3];
        lfm_float_re << float_re[2];
        lfm_float_re << float_re[1];
        lfm_float_re << float_re[0];

        lfm_float_im << float_im[3];
        lfm_float_im << float_im[2];
        lfm_float_im << float_im[1];
        lfm_float_im << float_im[0];
    }
    */

    ippsMul_32fc(data, data_shift, correl, fft_len);

    ippsFFTInv_CToC_32fc(correl, correl, pSpec, pDFTWorkBuf);
    

    for (size_t i = 0; i < fft_len; i++)
    {
        correl_re << i << '\t' << correl[i].re << std::endl;
    }

    for (size_t i = 0; i < fft_len; i++)
    {
        correl_im << i << '\t' << correl[i].im << std::endl;
    }

    ippsFree(data);
    ippsFree(pDFTSpec);
    ippsFree(pDFTInitBuf);
    ippsFree(pDFTWorkBuf);

	return 0;
}