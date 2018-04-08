#include <stdlib.h>
#include <stdio.h>
#include "bmp_reader.h"

#include "mpi.h"

void sequential_vertical_symmetry(int width, int_bmp_pixel_t (*tab)[width]){
        int i, j;
        for(i = 0; i < get_img_heigh(); i++)
        {
                for(j= 0; j < get_img_width()/2; j++)
                {
                        int_bmp_pixel_t temp = tab[i][j];
                        tab[i][j] = tab[i][get_img_width()-j-1];
                        tab[i][get_img_width()-j-1] = temp;
                }
        }
}

void sequential_horizontal_symmetry(int width, int_bmp_pixel_t (*tab)[width]){
        int i, j;
        for(i = 0; i < get_img_heigh()/2; i++)
        {
                for(j= 0; j < get_img_width(); j++)
                {
                        int_bmp_pixel_t temp = tab[i][j];
                        tab[i][j] = tab[get_img_heigh()-i-1][j];
                        tab[get_img_heigh()-i-1][j] = temp;
                }
        }
}

void parallel_vertical_symmetry(int width, int_bmp_pixel_t (*tab)[width], int nb_proc){
        int i, j;
        for(i = 0; i < get_img_heigh()/nb_proc; i++)
        {
                for(j= 0; j < get_img_width()/2; j++)
                {
                        int_bmp_pixel_t temp = tab[i][j];
                        tab[i][j] = tab[i][get_img_width()-j-1];
                        tab[i][get_img_width()-j-1] = temp;
                }
        }
}

int main(int argc, char* argv[])
{

        // int width = get_img_width_onemalloc("pingouin.bmp");
        // int_bmp_pixel_t (*tab)[width] = Lecture_image_onemalloc("pingouin.bmp");
        //
        // sequential_vertical_symmetry(width, tab);
        // sequential_horizontal_symmetry(width, tab);
        //
        // Ecriture_image_onemalloc(tab, "copie.bmp");
        // Liberation_image_lue_onemalloc(tab);
        // Liberation_finale();


        int rank, nbproc;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

        int width = get_img_width_onemalloc("pingouin.bmp");
        int height = get_img_heigh();

        int_bmp_pixel_t (*tab)[width] = Lecture_image_onemalloc("pingouin.bmp");
        int_bmp_pixel_t (*rcv)[width] = (int_bmp_pixel_t (*)[width])(int_bmp_pixel_t *) malloc((width)*(height)*sizeof(int_bmp_pixel_t));

        int nb_pixel_par_block = get_img_width()*get_img_heigh()/nbproc * 3;
        int sendcounts[nbproc];
        int displs[nbproc];

        for(int i = 0; i<nbproc; i++) {
                sendcounts[i] = nb_pixel_par_block;
                displs[i] = i * nb_pixel_par_block;
        }

        // divide the data among processes as described by sendcounts and displs
        MPI_Scatterv(tab, sendcounts, displs, MPI_INT, rcv, nb_pixel_par_block, MPI_INT, 0, MPI_COMM_WORLD);

        parallel_vertical_symmetry(width, rcv, nbproc);

        MPI_Gatherv(rcv, nb_pixel_par_block, MPI_INT, tab, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

        //transpose(width, height, tab);

        if(rank==0) {
                printf("%d %d\n", get_img_width(), get_img_heigh());
                Ecriture_image_onemalloc(tab, "copie.bmp");
        }


        Liberation_image_lue_onemalloc(tab);
        Liberation_image_lue_onemalloc(rcv);
        Liberation_finale();
        MPI_Finalize();

        return EXIT_SUCCESS;

}
