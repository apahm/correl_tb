#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <ipp.h>

const double Umax = 32767.0; //32767;
const double Umin = -32768.0;

int MakeMSequence(std::vector<int> poly, std::vector<int>& res)
{
    int k = poly.size() - 1;

    uint32_t N = std::pow(2, k) - 1;

    std::vector<std::vector<int>> matrix(k, std::vector<int>(k));

    matrix[0][k-1] = poly[0];
   
    for (int r = 1; r < k; r++)
    {
        matrix[r][r-1] = 1;
        matrix[r][k-1] = poly[r];
    }
    
    std::vector<int> X(k,0);
    X[0] = 1;
    
    std::vector<std::vector<int>> C(N, std::vector<int>(k));
    C[0][0] = 1;
    
    std::vector<int> Xtmp(k, 0);
    int j = 1;
    
    while (j < N)
    {

        for (int r = 0; r < k; r++)
        {
            Xtmp[r] = 0;
            for (int c = 0; c < k; c++)
                Xtmp[r] += matrix[r][c] * X[c];
        }
        for (int r = 0; r < k; r++)
        {
            X[r] = Xtmp[r] % 2;
            C[j][r] = X[r];
        }
        j++;
    }

    for (int r = 0; r < N; r++)
    {
        res.push_back(2 * C[r][0] - 1);
    }
    
    return 0;
}

int MakePolynom(std::vector<int> &poly, int max_degree)
{
    /*
        f(x) = x^13 + x^5 + x^2 + x + 1;
        N = 8191
    */
    if (max_degree == 13)
    {
        poly.push_back(1);
        poly.push_back(1);
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
    }
    /*
    f(x) = x^10 + x^3 + 1;
    N = 1023
    */
    else if (max_degree == 10)
    {
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
    }
    /*
    f(x) = x^8 + x^7 + x^2 + x + 1;
    N = 255
    */
    else if (max_degree == 8)
    {
        poly.push_back(1);
        poly.push_back(1);
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
        poly.push_back(1);
    }
    /*
    f(x) = x^7 + x + 1;
    N = 127
    */
    else if (max_degree == 7)
    {
        poly.push_back(1);
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
    }
    /*
    f(x) = x^6 + x + 1;
    N = 127
    */
    else if (max_degree == 6)
    {
        poly.push_back(1);
        poly.push_back(1);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(0);
        poly.push_back(1);
    }

    return 0;
}

int add_normal_distribution(Ipp32fc* data, double size, double SNR, double attenuator)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};

    for (size_t i = 0; i < size; i++)
    {
        data[i].re = attenuator * data[i].re;

        if (data[i].re > Umax)
            data[i].re = Umax;
        else if (data[i].re < -Umax)
            data[i].re = -Umax;

        data[i].im = attenuator * data[i].im;

        if (data[i].im > Umax)
            data[i].im = Umax;
        else if (data[i].im < -Umax)
            data[i].im = -Umax;
    }

    float sum_signal = 0.0;

    for (size_t i = 0; i < size; i++)
    {
        sum_signal += std::pow(data[i].re, 2);
    }

    double sigma_signal = sum_signal / size;

    double sigma_noise = sigma_signal / (std::pow(10, SNR / 10.0));

    std::normal_distribution<double> dist{0, std::sqrt(sigma_noise)};

    for (size_t i = 0; i < size; i++)
    {
        data[i].re = data[i].re + dist(gen);
        data[i].im = data[i].im + dist(gen);
    }

    float max_re = 0.0;
    float max_im = 0.0;

    for (size_t i = 1; i < size; i++)
    {
        if (data[i].re > max_re)
            max_re = data[i].re;
    }

    for (size_t i = 1; i < size; i++)
    {
        if (data[i].im > max_im)
            max_im = data[i].im;
    }

    for (size_t i = 0; i < size; i++)
    {
        data[i].re = Umax * (data[i].re / max_re);
        data[i].im = Umax * (data[i].im / max_im);
    }

    return 0;
}

int create_file_for_lab_work()
{
    std::ofstream lfm_int16_re("lfm_re_int16.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_int16_im("lfm_im_int16.dat", std::ios_base::out | std::ios_base::binary);

    std::ofstream lfm_int16_shift_re("lfm_int16_shift_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_int16_shift_im("lfm_int16_shift_im.dat", std::ios_base::out | std::ios_base::binary);

#ifdef LAB_WORK_PATH
    // Для лабораторного практикума - ЛЧМ - опорная функция
    std::ofstream lfm_re("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/lfm_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_im("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/lfm_im.dat", std::ios_base::out | std::ios_base::binary);

    // Для лабораторного практикума - сдвинутый ЛЧМ сигнал
    std::ofstream lfm_shift_re("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/lfm_shift_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_shift_im("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/lfm_shift_im.dat", std::ios_base::out | std::ios_base::binary);

    // Для лабораторного практикума - сдвинутый ЛЧМ сигнал с шумом
    std::ofstream lfm_shift_noise_re("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/lfm_shift_noise_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_shift_noise_im("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/lfm_shift_noise_im.dat", std::ios_base::out | std::ios_base::binary);

    // Для лабораторного практикума - результат корреляции
    std::ofstream correl_re("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/correl_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream correl_im("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/correl_im.dat", std::ios_base::out | std::ios_base::binary);
#else
    // Для лабораторного практикума - ЛЧМ - опорная функция
    std::ofstream lfm_re("lfm_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_im("lfm_im.dat", std::ios_base::out | std::ios_base::binary);

    // Для лабораторного практикума - сдвинутый ЛЧМ сигнал
    std::ofstream lfm_shift_re("lfm_shift_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_shift_im("lfm_shift_im.dat", std::ios_base::out | std::ios_base::binary);

    // Для лабораторного практикума - сдвинутый ЛЧМ сигнал с шумом
    std::ofstream lfm_shift_noise_re("lfm_shift_noise_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream lfm_shift_noise_im("lfm_shift_noise_im.dat", std::ios_base::out | std::ios_base::binary);

    // Для лабораторного практикума - результат корреляции
    std::ofstream correl_re("correl_re.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream correl_im("correl_im.dat", std::ios_base::out | std::ios_base::binary);
#endif
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

    // IPP_NODIV_BY_ANY
    ippsFFTGetSize_C_32fc(size, /*IPP_NODIV_BY_ANY*/IPP_NODIV_BY_ANY, ippAlgHintAccurate, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);

    pDFTSpec = ippsMalloc_8u(sizeDFTSpec);
    pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
    pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);

    ippsFFTInit_C_32fc(&pSpec, size,/*IPP_NODIV_BY_ANY*/ IPP_NODIV_BY_ANY, ippAlgHintAccurate, pDFTSpec, pDFTInitBuf);

    data = ippsMalloc_32fc(fft_len);
    data_shift = ippsMalloc_32fc(fft_len);

    ippsZero_32fc(data, fft_len);
    ippsZero_32fc(data_shift, fft_len);

    double lengthEtalonSignal = signalDuration * samplingFrequency;     // определение длинны эталонного сигнала
    for (int i = 0; i < lengthEtalonSignal; ++i)                        // формирование сигнала
    {
        double samplingT = (double)i / (double)samplingFrequency - (double)signalDuration / 2.0;    // Нахождение периода дискретизации

        // Формирование действительных составляющих
        double tempRe = round(Umax * cos(2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration));                       //TODO - заменить на cos 
        if (tempRe > Umax)
            tempRe = Umax;
        else if (tempRe < -Umax)
            tempRe = -Umax;

        // Формирование мнимых составляющих
        double tempIm = round(Umax * sin(2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration + fi * M_PI / 180.0));   //TODO - заменить на sin
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

    add_normal_distribution(data_shift, fft_len, -10, 1);

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

    for (size_t i = 0; i < fft_len; i++)
    {
        lfm_shift_noise_re << i << '\t' << data_shift[i].re << std::endl;
        lfm_shift_noise_im << i << '\t' << data_shift[i].im << std::endl;
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

int MakeLFMCenterFreq(double freqDev, double samplingFrequency, 
                    double center_freq, double signalDuration, 
                            double fi, std::vector<int16_t>& LFM, double sign)
{
    /*
    std::vector<int16_t> LFM;
    double freqDev = 1'000'000;
    double samplingFrequency = 100'000'000;
    double center_freq = 1'500'000;
    double signalDuration = 10 / 1000000.0;
    double fi = 0;
    double sign = 1;

    MakeLFMCenterFreq(freqDev, samplingFrequency, center_freq, signalDuration, fi, LFM, sign);
    */
    
    std::ofstream lfm_freq("C:/lab_fpga/Synopsis/data/third-party/lfm_freq.dat", std::ios_base::out | std::ios_base::binary);

    double lengthEtalonSignal = signalDuration * samplingFrequency;

    if (sign == 1)
    {
        for (int i = 0; i < lengthEtalonSignal; ++i)
        {
            double samplingT = (double)i / (double)samplingFrequency - (double)signalDuration / 2.0;

            double phase = 2.0 * M_PI * (double)center_freq * samplingT + (2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration);

            int16_t temp = round(Umax * cos(2.0 * M_PI * (double)center_freq * samplingT + (2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration)));
            if (temp > Umax)
                temp = Umax;
            else if (temp < -Umax)
                temp = -Umax;

            LFM.push_back(temp);
            lfm_freq << i << '\t' << temp << std::endl;
        }
    }
    else if (sign == -1)
    {
        for (int i = 0; i < lengthEtalonSignal; ++i)
        {
            double samplingT = (double)i / (double)samplingFrequency - (double)signalDuration / 2.0;

            double phase = 2.0 * M_PI * (double)center_freq * samplingT + (2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration);

            int16_t temp = round(Umax * cos(2.0 * M_PI * (double)center_freq * samplingT - (2.0 * M_PI * (double)freqDev * samplingT * samplingT / (double)signalDuration)));
            if (temp > Umax)
                temp = Umax;
            else if (temp < -Umax)
                temp = -Umax;

            LFM.push_back(temp);
            lfm_freq << i << '\t' << temp << std::endl;
        }
    }

    return 0;
}

enum class CSVState {
    UnquotedField,
    QuotedField,
    QuotedQuote
};

std::vector<std::string> readCSVRow(const std::string& row) {
    CSVState state = CSVState::UnquotedField;
    std::vector<std::string> fields{ "" };
    size_t i = 0; // index of the current field
    for (char c : row) {
        switch (state) {
        case CSVState::UnquotedField:
            switch (c) {
            case ',': // end of field
                fields.push_back(""); i++;
                break;
            case '"': state = CSVState::QuotedField;
                break;
            default:  fields[i].push_back(c);
                break;
            }
            break;
        case CSVState::QuotedField:
            switch (c) {
            case '"': state = CSVState::QuotedQuote;
                break;
            default:  fields[i].push_back(c);
                break;
            }
            break;
        case CSVState::QuotedQuote:
            switch (c) {
            case ',': // , after closing quote
                fields.push_back(""); i++;
                state = CSVState::UnquotedField;
                break;
            case '"': // "" -> "
                fields[i].push_back('"');
                state = CSVState::QuotedField;
                break;
            default:  // end of quote
                state = CSVState::UnquotedField;
                break;
            }
            break;
        }
    }
    return fields;
}

std::vector<std::vector<std::string>> readCSV(std::istream& in) {
    std::vector<std::vector<std::string>> table;
    std::string row;
    
    while (!in.eof()) {
        std::getline(in, row);
        if (in.bad() || in.fail()) {
            break;
        }
        auto fields = readCSVRow(row);
        table.push_back(fields);
    }
    return table;
}

int parse_ila_data()
{
    std::ifstream hardware_manager("C:/IP/iladata.csv", std::ios_base::in | std::ios_base::binary);

    std::vector<std::vector<std::string>> ila;

    if (hardware_manager.is_open())
    {
        ila = readCSV(hardware_manager);
    }

    std::vector<int16_t> adc_0_data;
    std::vector<int16_t> adc_1_data;

    for (size_t i = 1; i < ila.size(); i++)
    {
        uint8_t sample_adc_0[8];
        uint8_t sample_adc_1[8];

        std::string sample_str_adc_0;
        std::string sample_str_adc_1;

        for (size_t j = 0; j < 8; j++)
        {
            sample_str_adc_0.append(&ila[i][3][j * 2], 2);
            sample_str_adc_1.append(&ila[i][4][j * 2], 2);

            sample_adc_0[j] = std::stoi(sample_str_adc_0, nullptr, 16);
            sample_adc_1[j] = std::stoi(sample_str_adc_1, nullptr, 16);

            sample_str_adc_0.clear();
            sample_str_adc_1.clear();
        }

        int16_t sample_adc_0_0 = 0;
        int16_t sample_adc_0_1 = 0;
        int16_t sample_adc_0_2 = 0;
        int16_t sample_adc_0_3 = 0;

        std::memcpy((uint8_t*)(&sample_adc_0_3), &sample_adc_0[1], 1);
        std::memcpy((((uint8_t*)(&sample_adc_0_3)) + 1), &sample_adc_0[0], 1);
        std::memcpy((uint8_t*)(&sample_adc_0_2), &sample_adc_0[3], 1);
        std::memcpy((((uint8_t*)(&sample_adc_0_2)) + 1), &sample_adc_0[2], 1);
        std::memcpy((uint8_t*)(&sample_adc_0_1), &sample_adc_0[5], 1);
        std::memcpy((((uint8_t*)(&sample_adc_0_1)) + 1), &sample_adc_0[4], 1);
        std::memcpy((uint8_t*)(&sample_adc_0_0), &sample_adc_0[7], 1);
        std::memcpy((((uint8_t*)(&sample_adc_0_0)) + 1), &sample_adc_0[6], 1);

        int16_t sample_adc_1_0 = 0;
        int16_t sample_adc_1_1 = 0;
        int16_t sample_adc_1_2 = 0;
        int16_t sample_adc_1_3 = 0;

        std::memcpy((uint8_t*)(&sample_adc_1_3), &sample_adc_1[1], 1);
        std::memcpy((((uint8_t*)(&sample_adc_1_3)) + 1), &sample_adc_1[0], 1);
        std::memcpy((uint8_t*)(&sample_adc_1_2), &sample_adc_1[3], 1);
        std::memcpy((((uint8_t*)(&sample_adc_1_2)) + 1), &sample_adc_1[2], 1);
        std::memcpy((uint8_t*)(&sample_adc_1_1), &sample_adc_1[5], 1);
        std::memcpy((((uint8_t*)(&sample_adc_1_1)) + 1), &sample_adc_1[4], 1);
        std::memcpy((uint8_t*)(&sample_adc_1_0), &sample_adc_1[7], 1);
        std::memcpy((((uint8_t*)(&sample_adc_1_0)) + 1), &sample_adc_1[6], 1);

        adc_0_data.push_back(sample_adc_0_0);
        adc_0_data.push_back(sample_adc_0_1);
        adc_0_data.push_back(sample_adc_0_2);
        adc_0_data.push_back(sample_adc_0_3);

        adc_1_data.push_back(sample_adc_1_0);
        adc_1_data.push_back(sample_adc_1_1);
        adc_1_data.push_back(sample_adc_1_2);
        adc_1_data.push_back(sample_adc_1_3);
    }

    std::ofstream adc_0("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/adc_0.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream adc_1("C:/lab_fpga/lab_fpga/Synopsis/data/third-party/adc_1.dat", std::ios_base::out | std::ios_base::binary);

    for (size_t i = 0; i < ila.size(); i++)
    {
        adc_0 << i << '\t' << adc_0_data[i] << std::endl;
        adc_1 << i << '\t' << adc_1_data[i] << std::endl;
    }

    return 0;
}

int MakeMSequenceIQRecv(uint32_t max_degree, std::vector<float> &data_re, std::vector<float> &data_im)
{

    std::ofstream m_seq_re_recv("C:/lab_fpga/Synopsis/data/third-party/m_seq_re_recv.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream m_seq_im_recv("C:/lab_fpga/Synopsis/data/third-party/m_seq_im_recv.dat", std::ios_base::out | std::ios_base::binary);

    uint32_t N = std::pow(2, max_degree) - 1;

    std::vector<int> poly;
    std::vector<int> m_seq;

    MakePolynom(poly, max_degree);

    MakeMSequence(poly, m_seq);

    double freq_sample_rate = 400'000'000;
    double duration_code_interval = 10 / 1'000'000'000.0;

    uint32_t sample_per_coded_interval = duration_code_interval * freq_sample_rate;

    for (size_t i = 0; i < 380; i++)
    {
        data_re.push_back(0);
        data_im.push_back(0);
    }

    for (size_t i = 0; i < N; i++)
    {
        if (m_seq[i] == 1)
        {
            for (size_t j = 0; j < sample_per_coded_interval; j++)
            {
                data_re.push_back(Umax);
                data_im.push_back(0);
            }
        }
        else if (m_seq[i] == -1)
        {
            for (size_t j = 0; j < sample_per_coded_interval; j++)
            {
                data_re.push_back(Umin);
                data_im.push_back(0);
            }
        }
    }

    for (size_t i = 0; i < 2048 - 380 - N* sample_per_coded_interval; i++)
    {
        data_re.push_back(0);
        data_im.push_back(0);
    }

    for (size_t i = 0; i < data_re.size(); i++)
    {
        m_seq_re_recv << i << '\t' << data_re[i] << std::endl;
        m_seq_im_recv << i << '\t' << data_im[i] << std::endl;
    }

    return 0;
}

int MakeMSequenceIQSuppFunc(uint32_t max_degree, std::vector<float> &data_re, std::vector<float> &data_im)
{
    std::ofstream m_seq_re_sf("C:/lab_fpga/Synopsis/data/third-party/m_seq_re_sf.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream m_seq_im_sf("C:/lab_fpga/Synopsis/data/third-party/m_seq_im_sf.dat", std::ios_base::out | std::ios_base::binary);

    uint32_t N = std::pow(2, max_degree) - 1;

    std::vector<int> poly;
    std::vector<int> m_seq;

    MakePolynom(poly, max_degree);

    MakeMSequence(poly, m_seq);

    double freq_sample_rate = 400'000'000;
    double duration_code_interval = 10 / 1'000'000'000.0;

    uint32_t sample_per_coded_interval = duration_code_interval * freq_sample_rate;

    for (size_t i = 0; i < N; i++)
    {
        if (m_seq[i] == 1)
        {
            for (size_t j = 0; j < sample_per_coded_interval; j++)
            {
                data_re.push_back(Umax);
                data_im.push_back(0);
            }
        }
        else if (m_seq[i] == -1)
        {
            for (size_t j = 0; j < sample_per_coded_interval; j++)
            {
                data_re.push_back(Umin);
                data_im.push_back(0);
            }
        }
    }

    for (size_t i = 0; i < 2048 - N * sample_per_coded_interval; i++)
    {
        data_re.push_back(0);
        data_im.push_back(0);
    }

    for (size_t i = 0; i < data_re.size(); i++)
    {
        m_seq_re_sf << i << '\t' << data_re[i] << std::endl;
        m_seq_im_sf << i << '\t' << data_im[i] << std::endl;
    }

    return 0;
}

int main(int argc, char const *argv[])  
{


    // Для лабораторного практикума - результат корреляции
    std::ofstream correl_re("C:/lab_fpga/Synopsis/data/third-party/correl_re_m_seq.dat", std::ios_base::out | std::ios_base::binary);
    std::ofstream correl_im("C:/lab_fpga/Synopsis/data/third-party/correl_im_m_seq.dat", std::ios_base::out | std::ios_base::binary);

    uint32_t max_degree = 8;

    std::vector<float> data_re_recv;
    std::vector<float> data_im_recv;
    std::vector<float> data_re_sf;
    std::vector<float> data_im_sf;

    MakeMSequenceIQRecv(max_degree, data_re_recv, data_im_recv);

    MakeMSequenceIQSuppFunc(max_degree, data_re_sf, data_im_sf);

    int size = 11;

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

    // IPP_NODIV_BY_ANY
    ippsFFTGetSize_C_32fc(size, /*IPP_NODIV_BY_ANY*/IPP_NODIV_BY_ANY, ippAlgHintAccurate, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);

    pDFTSpec = ippsMalloc_8u(sizeDFTSpec);
    pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
    pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);

    ippsFFTInit_C_32fc(&pSpec, size,/*IPP_NODIV_BY_ANY*/ IPP_NODIV_BY_ANY, ippAlgHintAccurate, pDFTSpec, pDFTInitBuf);

    data = ippsMalloc_32fc(fft_len);
    data_shift = ippsMalloc_32fc(fft_len);

    ippsZero_32fc(data, fft_len);
    ippsZero_32fc(data_shift, fft_len);

    for (size_t i = 0; i < fft_len; i++)
    {
        data[i].re = data_re_sf[i];
        data[i].im = data_im_sf[i];
        
        data_shift[i].re = data_re_recv[i];
        data_shift[i].im = data_im_recv[i];
    }

    //add_normal_distribution(data_shift, fft_len, -10, 1);

    ippsFFTFwd_CToC_32fc(data, data, pSpec, pDFTWorkBuf);

    ippsFFTFwd_CToC_32fc(data_shift, data_shift, pSpec, pDFTWorkBuf);

    Ipp32fc* correl = nullptr;

    correl = ippsMalloc_32fc(fft_len);

    ippsZero_32fc(correl, fft_len);

    ippsConj_32fc(data, data, fft_len);

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