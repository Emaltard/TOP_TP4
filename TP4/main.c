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

void sequential_blur(int width, int_bmp_pixel_t (*tab)[width]){
        int_bmp_pixel_t (*temp)[width] = (int_bmp_pixel_t (*)[width])(int_bmp_pixel_t *) malloc((width)*(get_img_heigh())*sizeof(int_bmp_pixel_t));

        for(int i=0; i<get_img_heigh(); i++) {
                for(int j=0; j<width; j++) {
                        temp[i][j] = tab[i][j];
                }
        }

        for(int i=0; i<get_img_heigh(); i++) {
                for(int j=0; j<width; j++) {
                        int nbr_elem_for_blur = 0;
                        for(int k=i-1; k<=i+1; k++) {
                                for(int l=j-1; l<=j+1; l++) {
                                        if((k>=0) && (k<width) && (l>=0) && (l<get_img_heigh())) {
                                                if(!(k==i && l==j)) {
                                                        nbr_elem_for_blur++;
                                                        tab[i][j].Rouge += temp[k][l].Rouge;
                                                        tab[i][j].Bleu += temp[k][l].Bleu;
                                                        tab[i][j].Vert += temp[k][l].Vert;
                                                }
                                        }
                                }
                        }
                        tab[i][j].Rouge = tab[i][j].Rouge / nbr_elem_for_blur;
                        tab[i][j].Bleu = tab[i][j].Bleu / nbr_elem_for_blur;
                        tab[i][j].Vert = tab[i][j].Vert / nbr_elem_for_blur;
                }
        }

        Liberation_image_lue_onemalloc(temp);
}

void parallel_vertical_symmetry(int width, int_bmp_pixel_t (*tab)[width], int nbproc){

        int_bmp_pixel_t (*rcv)[width] = (int_bmp_pixel_t (*)[width])(int_bmp_pixel_t *) malloc((width)*(get_img_heigh())*sizeof(int_bmp_pixel_t));

        int nb_pixel_par_block = get_img_width()*get_img_heigh()/nbproc * 3;
        int sendcounts[nbproc];
        int displs[nbproc];

        for(int i = 0; i<nbproc; i++) {
                sendcounts[i] = nb_pixel_par_block;
                displs[i] = i * nb_pixel_par_block;
        }

        MPI_Scatterv(tab, sendcounts, displs, MPI_INT, rcv, nb_pixel_par_block, MPI_INT, 0, MPI_COMM_WORLD);

        int i, j;
        for(i = 0; i < get_img_heigh()/nbproc; i++)
        {
                for(j= 0; j < get_img_width()/2; j++)
                {
                        int_bmp_pixel_t temp = rcv[i][j];
                        rcv[i][j] = rcv[i][get_img_width()-j-1];
                        rcv[i][get_img_width()-j-1] = temp;
                }
        }

        MPI_Gatherv(rcv, nb_pixel_par_block, MPI_INT, tab, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

        Liberation_image_lue_onemalloc(rcv);
}

void parallel_horizontal_symmetry(int width, int_bmp_pixel_t (*tab)[width], int nbproc){
        int nb_pixel_par_block = get_img_width()/nbproc * 3;
        int sendcounts[nbproc];
        int displs[nbproc];

        for(int i = 0; i<nbproc; i++) {
                sendcounts[i] = nb_pixel_par_block;
                displs[i] = i * nb_pixel_par_block;
        }

        int_bmp_pixel_t (*rcv)[nb_pixel_par_block] = (int_bmp_pixel_t (*)[nb_pixel_par_block])(int_bmp_pixel_t *) malloc(nb_pixel_par_block*(get_img_heigh())*sizeof(int_bmp_pixel_t));


        for(int i=0; i < get_img_heigh(); i++) {
                MPI_Scatterv(tab[i], sendcounts, displs, MPI_INT, rcv[i], nb_pixel_par_block, MPI_INT, 0, MPI_COMM_WORLD);
        }

        for(int i = 0; i < get_img_heigh()/2; i++)
        {
                for(int j= 0; j < get_img_width()/nbproc; j++)
                {
                        int_bmp_pixel_t temp = rcv[i][j];
                        rcv[i][j] = rcv[get_img_heigh()-i-1][j];
                        rcv[get_img_heigh()-i-1][j] = temp;
                }
        }

        for(int i=0; i < get_img_heigh(); i++) {
                MPI_Gatherv(rcv[i], nb_pixel_par_block, MPI_INT, tab[i], sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
        }

        Liberation_image_lue_onemalloc(rcv);
}

int main(int argc, char* argv[])
{
        int rank, nbproc;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

        int width = get_img_width_onemalloc("pingouin.bmp");

        int_bmp_pixel_t (*tab)[width] = Lecture_image_onemalloc("pingouin.bmp");

        parallel_vertical_symmetry(width, tab, nbproc);

        parallel_horizontal_symmetry(width, tab, nbproc);


        if(rank==0) {
                int_bmp_pixel_t (*tab2)[width] = Lecture_image_onemalloc("pingouin.bmp");
                sequential_blur(width, tab2);
                printf("%d %d\n", get_img_width(), get_img_heigh());
                Ecriture_image_onemalloc(tab, "copie.bmp");
                Ecriture_image_onemalloc(tab2, "copie2.bmp");
                Liberation_image_lue_onemalloc(tab2);
        }

        Liberation_image_lue_onemalloc(tab);
        Liberation_finale();
        MPI_Finalize();

        return EXIT_SUCCESS;

}
