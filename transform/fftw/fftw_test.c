//
// Created by wcy on 2023/2/18.
//

//#include <sys/time.h>
#include <stdlib.h>
#include <fftw3.h>
#include <time.h>
#include <math.h>


int main() {
    fftw_complex *in, *out, *h;
    fftw_complex *buffer_in, *buffer_out;
    int N = 1 << 15;
    fftw_complex buff;

    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    h = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    buffer_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    buffer_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = 1 + i;
        in[i][1] = 2 + i;
    }
    for (int i = N / 2; i < N; i++) {
//        in[i] = 0;
//        h[i] = 0;
    }

    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan q = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan r = fftw_plan_dft_1d(N, buffer_out, buffer_in, FFTW_FORWARD, FFTW_ESTIMATE);

//    fftw_export_wisdom_to_filename("wisdom.txt");

    // start timer
//    clock_t start = clock(), diff;
//    struct timeval start, end;
//    gettimeofday(&start, NULL);

    fftw_execute(p);
    fftw_execute(q);

    // stop timer
//    diff = clock() - start;
//    gettimeofday(&end, NULL);

    fftw_execute(q);
    for (int i = 0; i < N; ++i) {
        buff[0] = buffer_out[i][0] * buffer_out[i][0] - buffer_out[i][1] * buffer_out[i][1];
        buff[1] = buffer_in[i][0] * buffer_out[i][1] + buffer_out[i][0] * buffer_in[i][1];
        buffer_in[i][0] = buff[0];
        buffer_in[i][1] = buff[1];
    }

    fftw_execute(r);
    for (int i = 0; i < N; ++i) {
//        out[i] /= N;
    }

    // get end time
//    clock_t end = clock();


    for (int i = 0; i < 10; i++) {
        printf("%f + %fi\n", in[i][0] / N, in[i][1] / N);
    }

//    printf("Time taken %lf ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    fftw_destroy_plan(p);
    fftw_destroy_plan(q);
    fftw_destroy_plan(r);
    fftw_free(in);
    fftw_free(out);
    fftw_free(h);
    fftw_free(buffer_in);
    fftw_free(buffer_out);

}


void fftw_main_r() {
    int N = 1 << 15;
    double *in;
    fftw_complex *out;
    fftw_plan p, q;

    in = (double *) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);

    //random
    srand(time(NULL));

    for (int i = 0; i < N; i++) {
        in[i] = 1 + i + rand() % 10;
    }

    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    q = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);


    fftw_execute(p);
    fftw_execute(q);
    // stop timer




    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

void fftw_sin() {
//    int N = 1 << 20;
//    double *in;
//    fftw_complex *out;
//    fftw_plan p;
//
//    in = (double *) fftw_malloc(sizeof(double) * N);
//    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
//
//    for (int i = 0; i < N; i++) {
//        in[i] = sin(2 * M_PI * 1000 * i / N);
//    }

}
