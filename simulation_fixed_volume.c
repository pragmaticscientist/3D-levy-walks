#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> // Per strcmp e per gestire le stringhe
#include <sys/stat.h>
#include <errno.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "func.h"
#include "helpers.h"

double get_diameter_from_measure_and_shape(double measure, const char* shape) {
    if (measure <= 0) {
        fprintf(stderr, "Error: The measure (volume/area/length) must be positive.\n");
        return -1.0;
    }

    if (strcmp(shape, "Ball") == 0) {
        // Measure is 3D Volume for a Sphere
        // V = (1/6) * PI * D^3
        // D = cbrt(6 * V / PI)
        return cbrt(6.0 * measure / M_PI);
    } else if (strcmp(shape, "Disk") == 0) {
        // Measure is 2D Area for a Circle (Disk)
        // A = (1/4) * PI * D^2
        // D = sqrt(4 * A / PI)
        return sqrt(4.0 * measure / M_PI);
    } else if (strcmp(shape, "Line") == 0) {
        // Measure is 1D Length for a Line
        // L = D
        // V = pi * D + 4/3 pi
        return measure;
    } else {
        fprintf(stderr, "Error: Unknown shape '%s'. Supported shapes are 'Ball', 'Disk', 'Line'.\n", shape);
        return -1.0; // Indicate error
    }
}

static int ensure_directory_exists(const char *path) {
    if (!path || !*path) {
        return 0;
    }

    char *copy = strdup(path);
    if (!copy) {
        return -1;
    }

    for (char *p = copy + 1; *p; ++p) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(copy, 0777) != 0 && errno != EEXIST) {
                free(copy);
                return -1;
            }
            *p = '/';
        }
    }

    if (mkdir(copy, 0777) != 0 && errno != EEXIST) {
        free(copy);
        return -1;
    }

    free(copy);
    return 0;
}

static char *build_thread_output_path(const char *save_dir, const char *file_name, int thread_index) {
    if (!save_dir || !file_name) {
        return NULL;
    }

    const char *dot = strrchr(file_name, '.');
    size_t base_len = strlen(file_name);
    size_t dir_len = strlen(save_dir);
    size_t suffix_len = 32;
    size_t ext_len = dot ? strlen(dot) : 0;

    char *out = malloc(dir_len + base_len + suffix_len + 4);
    if (!out) {
        return NULL;
    }

    if (save_dir[dir_len - 1] == '/') {
        if (dot) {
            size_t stem_len = (size_t)(dot - file_name);
            snprintf(out, dir_len + base_len + suffix_len + 4, "%s%.*s_thread%d%s", save_dir, (int)stem_len, file_name, thread_index, dot);
        } else {
            snprintf(out, dir_len + base_len + suffix_len + 4, "%s%s_thread%d", save_dir, file_name, thread_index);
        }
    } else {
        if (dot) {
            size_t stem_len = (size_t)(dot - file_name);
            snprintf(out, dir_len + base_len + suffix_len + 4, "%s/%.*s_thread%d%s", save_dir, (int)stem_len, file_name, thread_index, dot);
        } else {
            snprintf(out, dir_len + base_len + suffix_len + 4, "%s/%s_thread%d", save_dir, file_name, thread_index);
        }
    }

    (void)ext_len;
    return out;
}

static int write_combined_output(const char *output_path, const char *header, FILE **thread_files, int num_threads) {
    FILE *combined = fopen(output_path, "w");
    if (!combined) {
        perror("Error opening combined output file");
        return -1;
    }

    fputs(header, combined);
    for (int tid = 0; tid < num_threads; ++tid) {
        if (!thread_files[tid]) {
            continue;
        }
        if (fflush(thread_files[tid]) != 0) {
            perror("Error flushing thread output file");
            fclose(combined);
            return -1;
        }
        rewind(thread_files[tid]);
        char line[4096];
        while (fgets(line, sizeof(line), thread_files[tid])) {
            if (strcmp(line, header) == 0) {
                continue;
            }
            fputs(line, combined);
        }
    }

    fclose(combined);
    return 0;
}

int main(int argc, char *argv[]) {
    double fixed_target_dist = 0;  // Default: target on origin


    if (argc != 2) {
        fprintf(stderr, "Error! Add a configuration file.\n");
        return EXIT_FAILURE;
    }

    const char* conf_path = argv[1];

    // load configuration
    Config cfg = {0};
    load_config(conf_path, &cfg);

    printf("## Starting Levy Search Simulations on 3D Torus\n");

    if (ensure_directory_exists(cfg.save_directory) != 0) {
        fprintf(stderr, "Error creating output directory: %s\n", cfg.save_directory);
        return EXIT_FAILURE;
    }

    char *output_filename = build_output_path(cfg.save_directory, cfg.file_name);
    const char *csv_header = "n_walkers,n_volume,mu,lmax,D,TargetShape,DiskDiameter,n_targets,fixed_target_dist,detection_time,probability\n";

    int num_threads = cfg.num_threads > 0 ? cfg.num_threads : 1;
#ifdef _OPENMP
    if (num_threads > 1) {
        omp_set_num_threads(num_threads);
    }
#endif

    FILE **thread_files = NULL;
    char **thread_paths = NULL;
    if (num_threads > 1) {
        thread_files = calloc((size_t)num_threads, sizeof(FILE *));
        thread_paths = calloc((size_t)num_threads, sizeof(char *));
        if (!thread_files || !thread_paths) {
            fprintf(stderr, "Error allocating thread output buffers.\n");
            free(output_filename);
            free_config(&cfg);
            return EXIT_FAILURE;
        }

        for (int tid = 0; tid < num_threads; ++tid) {
            thread_paths[tid] = build_thread_output_path(cfg.save_directory, cfg.file_name, tid);
            if (!thread_paths[tid]) {
                fprintf(stderr, "Error building output path for thread %d.\n", tid);
                for (int i = 0; i < tid; ++i) {
                    free(thread_paths[i]);
                    fclose(thread_files[i]);
                }
                free(thread_paths);
                free(thread_files);
                free(output_filename);
                free_config(&cfg);
                return EXIT_FAILURE;
            }

            thread_files[tid] = fopen(thread_paths[tid], "w+");
            if (!thread_files[tid]) {
                perror("Error opening thread output file");
                for (int i = 0; i < tid; ++i) {
                    free(thread_paths[i]);
                    fclose(thread_files[i]);
                }
                free(thread_paths[tid]);
                free(thread_paths);
                free(thread_files);
                free(output_filename);
                free_config(&cfg);
                return EXIT_FAILURE;
            }

            fputs(csv_header, thread_files[tid]);
        }
        printf("Using %d OpenMP threads; per-thread CSV files will be written under %s\n", num_threads, cfg.save_directory);
    } else {
        FILE *output_file = fopen(output_filename, "w");
        if (output_file == NULL) {
            perror("Error opening output file");
            free(output_filename);
            free_config(&cfg);
            return EXIT_FAILURE;
        }
        fputs(csv_header, output_file);
        fclose(output_file);
    }

    long start_total_time = time(NULL);

    long total_inner_iterations = (long)cfg.len_range_side * cfg.len_list_shapes * cfg.len_range_disk_diameter *
                                   cfg.len_range_ntargets * cfg.len_range_nwalkers * cfg.len_rangemu_LevyDistrib * cfg.len_range_probability;
    
    printf("Warning: this code considers the volume of the extended target (target + boundary) for the calculation of D.\n");

    if (num_threads > 1) {
#ifdef _OPENMP
#pragma omp parallel num_threads(num_threads)
        {
            int tid = omp_get_thread_num();
            FILE *thread_output = thread_files[tid];
            #pragma omp for schedule(static)
            for (int current_trial = 0; current_trial < cfg.num_runs; ++current_trial) {
                long pbar_counter = 0;
                int temp1;
                double temp2;

                for (int i_side = 0; i_side < cfg.len_range_side; ++i_side) {
                    int side = cfg.range_side[i_side];
                    double n_volume = pow(side, 3);
                    int lmax_current = side / 2;
                    for (int i_mu = 0; i_mu < cfg.len_rangemu_LevyDistrib; ++i_mu) {
                        double mu = cfg.rangemu_LevyDistrib[i_mu];
                        double a = get_normalization_constant(mu, lmax_current);
                        for (int i_disk_diameter = 0; i_disk_diameter < cfg.len_range_disk_diameter; ++i_disk_diameter) {
                            double disk_diameter = cfg.range_disk_diameter[i_disk_diameter];
                            for (int i_ntargets = 0; i_ntargets < cfg.len_range_ntargets; ++i_ntargets) {
                                int n_targets = cfg.range_ntargets[i_ntargets];
                                for (int i_nwalkers = 0; i_nwalkers < cfg.len_range_nwalkers; ++i_nwalkers) {
                                    int n_walkers = cfg.range_nwalkers[i_nwalkers];
                                    for (int i_shape = 0; i_shape < cfg.len_list_shapes; ++i_shape) {
                                        const char *TargetShape = cfg.list_shapes[i_shape];
                                        double D;
                                        if (strcmp(TargetShape, cfg.reference_shape) == 0) {
                                            D = disk_diameter;
                                        } else if (strcmp(TargetShape, "Ball_no_boundary") == 0) {
                                            double target_volume = 4 * M_PI + 2 * disk_diameter * M_PI;
                                            D = get_diameter_from_measure_and_shape(target_volume, "Ball");
                                        } else {
                                            fprintf(stderr, "Error: Unsupported shape '%s'.\n", TargetShape);
                                            continue;
                                        }
                                        for (int i_prob = 0; i_prob < cfg.len_range_probability; ++i_prob) {
                                            double p = cfg.range_probability[i_prob];
                                            Result result = LevySearch3D_MultiWalker(n_walkers, "nest", n_volume, mu, lmax_current,
                                                                                         D, TargetShape, n_targets, fixed_target_dist, p, a, 2, 0, 0, 0.0);
                                            double detection_time = result.detection_time;
                                            fprintf(thread_output, "%d,%.0f,%.1f,%d,%.2f,%s,%.0f,%d,%.1f,%.2f,%.2f\n",
                                                    n_walkers, n_volume, mu, lmax_current, D, TargetShape, disk_diameter, n_targets, fixed_target_dist, detection_time, p);
                                            pbar_counter++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
#else
        (void)num_threads;
#endif
    } else {
        FILE *output_file = fopen(output_filename, "w");
        if (output_file == NULL) {
            perror("Error opening output file");
            free(output_filename);
            free_config(&cfg);
            return EXIT_FAILURE;
        }
        fputs(csv_header, output_file);
        for (int current_trial = 0; current_trial < cfg.num_runs; ++current_trial) {
            printf("\n--- Overall Trial: %d/%d ---\n", current_trial + 1, cfg.num_runs);
            long elapsed_seconds = time(NULL) - start_total_time;
            printf("    %d trials completed, total elapsed time: %.2f minutes\n", current_trial + 1, (double)elapsed_seconds / 60.0);
            if (current_trial > 0) {
                double estimated_remaining_time = ((double)elapsed_seconds / (current_trial + 1)) * (cfg.num_runs - (current_trial + 1));
                printf("    Estimated time remaining: %.2f minutes\n", estimated_remaining_time / 60.0);
            }

            long pbar_counter = 0;
            int temp1;
            double temp2;

            for (int i_side = 0; i_side < cfg.len_range_side; ++i_side) {
                int side = cfg.range_side[i_side];
                double n_volume = pow(side, 3);
                int lmax_current = side / 2;
                for (int i_mu = 0; i_mu < cfg.len_rangemu_LevyDistrib; ++i_mu) {
                    double mu = cfg.rangemu_LevyDistrib[i_mu];
                    double a = get_normalization_constant(mu, lmax_current);
                    for (int i_disk_diameter = 0; i_disk_diameter < cfg.len_range_disk_diameter; ++i_disk_diameter) {
                        double disk_diameter = cfg.range_disk_diameter[i_disk_diameter];
                        for (int i_ntargets = 0; i_ntargets < cfg.len_range_ntargets; ++i_ntargets) {
                            int n_targets = cfg.range_ntargets[i_ntargets];
                            for (int i_nwalkers = 0; i_nwalkers < cfg.len_range_nwalkers; ++i_nwalkers) {
                                int n_walkers = cfg.range_nwalkers[i_nwalkers];
                                for (int i_shape = 0; i_shape < cfg.len_list_shapes; ++i_shape) {
                                    const char *TargetShape = cfg.list_shapes[i_shape];
                                    double D;
                                    if (strcmp(TargetShape, cfg.reference_shape) == 0) {
                                        D = disk_diameter;
                                    } else if (strcmp(TargetShape, "Ball_no_boundary") == 0) {
                                        double target_volume = 4 * M_PI + 2 * disk_diameter * M_PI;
                                        D = get_diameter_from_measure_and_shape(target_volume, "Ball");
                                    } else {
                                        fprintf(stderr, "Error: Unsupported shape '%s'.\n", TargetShape);
                                        continue;
                                    }
                                    for (int i_prob = 0; i_prob < cfg.len_range_probability; ++i_prob) {
                                        double p = cfg.range_probability[i_prob];
                                        Result result = LevySearch3D_MultiWalker(n_walkers, "nest", n_volume, mu, lmax_current,
                                                                                     D, TargetShape, n_targets, fixed_target_dist, p, a, 2, 0, 0, 0.0);
                                        double detection_time = result.detection_time;
                                        fprintf(output_file, "%d,%.0f,%.1f,%d,%.2f,%s,%.0f,%d,%.1f,%.2f,%.2f\n",
                                                n_walkers, n_volume, mu, lmax_current, D, TargetShape, disk_diameter, n_targets, fixed_target_dist, detection_time, p);
                                        pbar_counter++;
                                        if (total_inner_iterations > 0 && pbar_counter % (total_inner_iterations / 100 + 1) == 0) {
                                            printf("\rTrial %d/%d Progress: %.2f%%", current_trial + 1, cfg.num_runs,
                                                   (double)pbar_counter * 100.0 / total_inner_iterations);
                                            fflush(stdout);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            printf("\rTrial %d/%d Progress: 100.00%%\n", current_trial + 1, cfg.num_runs);
        }
        fclose(output_file);
    }

    long end_total_time = time(NULL);
    printf("\nTotal simulation time: %.2f minutes\n", (double)(end_total_time - start_total_time) / 60.0);

    if (num_threads > 1) {
        if (write_combined_output(output_filename, csv_header, thread_files, num_threads) != 0) {
            fprintf(stderr, "Error writing combined output file.\n");
        } else {
            printf("Combined results written to %s\n", output_filename);
        }
    }

    for (int tid = 0; tid < num_threads; ++tid) {
        if (thread_files && thread_files[tid]) {
            fclose(thread_files[tid]);
        }
        free(thread_paths ? thread_paths[tid] : NULL);
    }
    free(thread_paths);
    free(thread_files);
    free_config(&cfg); 
    printf("Saving results to %s\n", output_filename);
    free(output_filename);
    printf("Save complete.\n");

    return 0;
}
