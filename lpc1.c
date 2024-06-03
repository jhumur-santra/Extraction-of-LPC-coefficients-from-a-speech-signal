#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// pre-emphasis
void pre_emphasis(float *signal, int length, float a) {
    for (int i = length - 1; i > 0; i--) {
        signal[i] -= a * signal[i - 1];
    }
}

// frame blocking
float** frame_blocking(float *signal, int M, int N, int L) {
    float **frames = (float **)malloc(L * sizeof(float *));
    if (frames == NULL) {
        printf("Memory allocation error\n");
        return NULL;
    }

    for (int l = 0; l < L; l++) {
        frames[l] = (float *)malloc(N * sizeof(float));
        if (frames[l] == NULL) {
            printf("Memory allocation error\n");
            for (int i = 0; i < l; i++) {
                free(frames[i]);
            }
            free(frames);
            return NULL;
        }

        for (int n = 0; n < N; n++) {
            frames[l][n] = signal[(M * l) + n];
        }
    }

    return frames;
}

// windowing (hamming window)
void windowing(float *frames, int len) {
    for (int i = 0; i < len; i++) {
        frames[i] *= 0.54 - (0.46 * cos(2 * M_PI * i / (len - 1)));
    }
}

// autocorrelation analysis
float *autocorrelation(float *frames, int p, int N) {
    float *r = (float *)calloc(p + 1, sizeof(float));
    if (r == NULL) {
        printf("Memory allocation error\n");
        return NULL;
    }

    for (int m = 0; m <= p; m++) {
        r[m] = 0;
        for (int n = 0; n <= N - 1 - m; n++) {
            r[m] += frames[n] * frames[n + m];
        }
    }
    return r;
}

// LPC analysis
void lpc(float *r, float *a, float *k, int p) {
    float *E = (float *)malloc((p + 1) * sizeof(float));
    float *alpha_prev = (float *)malloc((p + 1) * sizeof(float));
    float *alpha_curr = (float *)malloc((p + 1) * sizeof(float));

    E[0] = r[0];

    for (int i = 1; i <= p; i++) {
        float sum = 0.0;
        for (int j = 1; j < i; j++) {
            sum += alpha_prev[j] * r[i - j];
        }

        k[i] = (r[i] - sum) / E[i - 1];
        alpha_curr[i] = k[i];

        for (int j = 1; j < i; j++) {
            alpha_curr[j] = alpha_prev[j] - k[i] * alpha_prev[i - j];
        }

        E[i] = (1 - k[i] * k[i]) * E[i - 1];

        for (int j = 1; j <= i; j++) {
            alpha_prev[j] = alpha_curr[j];
        }
    }

    for (int i = 1; i <= p; i++) {
        a[i] = alpha_curr[i];
    }
}

// main function
int main() {
    FILE *file;
    float *signal = NULL;
    int num;
    int count = 0;

    file = fopen("234101011_a_1.txt", "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    while (fscanf(file, "%d", &num) == 1) {
        float *temp = realloc(signal, (count + 1) * sizeof(float));
        if (temp == NULL) {
            printf("Memory allocation error\n");
            free(signal);
            fclose(file);
            return 1;
        }
        signal = temp;

        signal[count] = num;
        count++;
    }

    fclose(file);

    float a = 0.95; // Pre-emphasis factor
    pre_emphasis(signal, count, a);

    int M = 100;
    int N = 653;
    int L = count / N;
    float **frames = frame_blocking(signal, M, N, L);

    for (int l = 0; l < L; l++) {
        windowing(frames[l], N);
    }

    int p = 10;
    for (int l = 0; l < L; l++) {
        float *r = autocorrelation(frames[l], p, N);
        float *a = (float *)malloc((p + 1) * sizeof(float));
        float *k = (float *)malloc((p + 1) * sizeof(float));

        lpc(r, a, k, p);

        printf("Frame %d:\n", l);
        printf("LPC Coefficients:\n");
        for (int i = 1; i <= p; i++) {
            printf("a[%d] = %f\n", i, a[i]);
        }

        printf("\nPARCOR Coefficients:\n");
        for (int i = 1; i <= p; i++) {
            printf("k[%d] = %f\n", i, k[i]);
        }
    }

    return 0;
}

