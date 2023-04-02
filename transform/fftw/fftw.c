//
// Created by wcy on 2023/2/18.
//

#include <time.h>
#include <sys/time.h>
#include "fftw.h"


void fftw_main() {
    fftw_complex *in, *out, *h;
    fftw_complex *buffer_in, *buffer_out;
    int N = 1 << 14;
    fftw_complex buff;

    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    h = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    buffer_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    buffer_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N / 2; i++) {
        in[i][0] = 1 + i;
//        h[i] = 1 + i;
    }
    for (int i = N / 2; i < N; i++) {
//        in[i] = 0;
//        h[i] = 0;
    }

    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan q = fftw_plan_dft_1d(N, buffer_in, buffer_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan r = fftw_plan_dft_1d(N, buffer_out, buffer_in, FFTW_FORWARD, FFTW_ESTIMATE);

//    fftw_export_wisdom_to_filename("wisdom.txt");
    // get start time
    struct timeval start, end;
    gettimeofday(&start, NULL);

    fftw_execute(p);
    fftw_execute(q);

    // get end time
    gettimeofday(&end, NULL);
    printf("fftw time: %f ms\n",
           (double) (end.tv_sec - start.tv_sec) * 1000 + (double) (end.tv_usec - start.tv_usec) / 1000);

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
        printf("%f\n", out[i]);
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(q);
    fftw_destroy_plan(r);
    fftw_free(in);
    fftw_free(out);
    fftw_free(h);
    fftw_free(buffer_in);
    fftw_free(buffer_out);

}
