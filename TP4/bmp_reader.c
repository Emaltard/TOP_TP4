#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "bmp_internal.h"
#include "bmp_reader.h"

#define ONEALLOC "ONEALLOC"

static bmp_header_infos_t *header = NULL;

static void
dump_header_bmp(header_bmp_t* img)
{
        fprintf(stdout, "TailleFichier: %d\n", img->TailleFichier);
        fprintf(stdout, "Reserves: %d\n", img->Reserves);
        fprintf(stdout, "PositionImage: %d\n", img->PositionImage);
        fprintf(stdout, "TailleEntete: %d\n", img->TailleEntete);
        fprintf(stdout, "Width: %d\n", img->Width);
        fprintf(stdout, "Height: %d\n", img->Height);
        fprintf(stdout, "NombrePlan: %d\n", img->NombrePlan);
        fprintf(stdout, "NombreBit: %d\n", img->NombreBit);
        fprintf(stdout, "TypeCompres: %d\n", img->TypeCompres);
        fprintf(stdout, "TailleImage: %d\n", img->TailleImage);
        fprintf(stdout, "ResolutionH: %d\n", img->ResolutionH);
        fprintf(stdout, "ResolutionV: %d\n", img->ResolutionV);
        fprintf(stdout, "NombreCouleur: %d\n", img->NombreCouleur);
        fprintf(stdout, "CouleurImport: %d\n", img->CouleurImport);
        return;
}

static void
readBMP_header(FILE *img)
{
        char FileType[3];               // Obligatoirement Â« BM Â» pour du BMP

        /* Allocate bmp_header */
        if (header == NULL)
                header = (bmp_header_infos_t *) malloc(sizeof(bmp_header_infos_t));
        assert(header != NULL);


        /* Check img filetype */
        if(fread(FileType, sizeof(char), 2, img)) {
                FileType[2] = '\0';
                assert(strcmp(FileType, "BM") == 0);
        }


        /* Check img header */
        if(!fread(&(header->raw), sizeof(header_bmp_t), 1, img)) {
                abort();
        }

        /* Check img header size and first pixel position */
        assert( (sizeof(header_bmp_t) + 2) == header->raw.PositionImage);

        /* Extract image height and width from header */
        header->width   = header->raw.Width;
        header->height  = header->raw.Height;

}

static int_bmp_pixel_t *
readBMP_pixels_onemalloc(FILE *img)
{
        int i, j, ret;
        int_bmp_pixel_t * PixelsTab = NULL;
        char_bmp_pixel_t tmp;
        assert(header != NULL);
        PixelsTab = (int_bmp_pixel_t *) malloc((header->width)*(header->height)*sizeof(int_bmp_pixel_t));
        assert(PixelsTab != NULL);

        ret = 0;
        for(i = 0; i < header->height; i++)
        {
                for(j = 0; j < header->width; j++)
                {
                        ret += fread(&tmp, sizeof(char_bmp_pixel_t), 1, img);
                        PixelsTab[i*(header->width)+j].Rouge = (tmp.Rouge);
                        //assert(PixelsTab_p[i][j].Rouge >= 0 && PixelsTab_p[i][j].Rouge < 256);
                        PixelsTab[i*(header->width)+j].Bleu = (int) (tmp.Bleu);
                        //assert(PixelsTab_p[i][j].Bleu >= 0 && PixelsTab_p[i][j].Bleu < 256);
                        PixelsTab[i*header->width+j].Vert = (int) (tmp.Vert);
                        //assert(PixelsTab[i][j].Vert >= 0 && PixelsTab[i][j].Vert < 256);
                }
        }

        assert(ret == header->height*header->width);

        return PixelsTab;
}

static int_bmp_pixel_t **
readBMP_pixels(FILE *img)
{
        int i, j, ret;
        int_bmp_pixel_t ** PixelsTab = NULL;
        char_bmp_pixel_t tmp;
        assert(header != NULL);
        PixelsTab = (int_bmp_pixel_t **) malloc(header->height*sizeof(int_bmp_pixel_t *));
        assert(PixelsTab != NULL);

        for(i = 0; i < header->height; i++)
        {
                PixelsTab[i] = (int_bmp_pixel_t *) malloc(header->width*sizeof(int_bmp_pixel_t));
                assert(PixelsTab[i] != NULL);
        }

        ret = 0;
        for(i = 0; i < header->height; i++)
        {
                for(j = 0; j < header->width; j++)
                {
                        ret += fread(&tmp, sizeof(char_bmp_pixel_t), 1, img);
                        PixelsTab[i][j].Rouge = (tmp.Rouge);
                        //assert(PixelsTab_p[i][j].Rouge >= 0 && PixelsTab_p[i][j].Rouge < 256);
                        PixelsTab[i][j].Bleu = (int) (tmp.Bleu);
                        //assert(PixelsTab_p[i][j].Bleu >= 0 && PixelsTab_p[i][j].Bleu < 256);
                        PixelsTab[i][j].Vert = (int) (tmp.Vert);
                        //assert(PixelsTab[i][j].Vert >= 0 && PixelsTab[i][j].Vert < 256);
                }
        }

        assert(ret == header->height*header->width);

        return (int_bmp_pixel_t **)PixelsTab;
}

static void
writeBMP_header(FILE* img)
{
        char *FileType = "BM";

        /* write BMP filetype */
        fwrite(FileType, sizeof(char), 2, img);

        /* copy header info */
        fwrite(&(header->raw), sizeof(header_bmp_t), 1, img);

        return;
}

static void
writeBMP_pixels_onemalloc( int_bmp_pixel_t (*PixelsTab)[header->width], FILE* img)
{
        int i, j, ret;
        char_bmp_pixel_t tmp;

        ret = 0;
        for(i = 0; i < header->height; i++)
        {
                for(j = 0; j < header->width; j++)
                {
                        tmp.Rouge = (char) (PixelsTab[i][j].Rouge);
                        //assert( PixelsTab[i][j].Rouge >= 0 &&  PixelsTab[i][j].Rouge < 256);
                        tmp.Bleu = (char) (PixelsTab[i][j].Bleu);
                        //assert( PixelsTab[i][j].Bleu >= 0 &&  PixelsTab[i][j].Bleu < 256);
                        tmp.Vert = (char) (PixelsTab[i][j].Vert);
                        //assert( PixelsTab[i][j].Vert >= 0 &&  PixelsTab[i][j].Vert < 256);
                        ret += fwrite(&tmp, sizeof(char_bmp_pixel_t), 1, img);
                }
        }
        assert(ret == header->height*header->width);
}
static void
writeBMP_pixels( int_bmp_pixel_t ** PixelsTab, FILE* img)
{
        int i, j, ret;
        char_bmp_pixel_t tmp;

        ret = 0;
        for(i = 0; i < header->height; i++)
        {
                for(j = 0; j < header->width; j++)
                {
                        tmp.Rouge = (char) (PixelsTab[i][j].Rouge);
                        //assert( PixelsTab[i][j].Rouge >= 0 &&  PixelsTab[i][j].Rouge < 256);
                        tmp.Bleu = (char) (PixelsTab[i][j].Bleu);
                        //assert( PixelsTab[i][j].Bleu >= 0 &&  PixelsTab[i][j].Bleu < 256);
                        tmp.Vert = (char) (PixelsTab[i][j].Vert);
                        //assert( PixelsTab[i][j].Vert >= 0 &&  PixelsTab[i][j].Vert < 256);
                        ret += fwrite(&tmp, sizeof(char_bmp_pixel_t), 1, img);
                }
        }
        assert(ret == header->height*header->width);
}

static int
writeBMP_img_onemalloc(int_bmp_pixel_t (*pxl)[], char* filename)
{
        FILE *img = NULL;

        /* Openfile */
        img = fopen(filename, "w");
        assert(img != NULL);

        writeBMP_header(img);
        writeBMP_pixels_onemalloc(pxl, img);

        /* close file */
        fclose(img);
        return 0;
}
static int
writeBMP_img(int_bmp_pixel_t ** pxl, char* filename)
{
        FILE *img = NULL;

        /* Openfile */
        img = fopen(filename, "w");
        assert(img != NULL);

        writeBMP_header(img);
        writeBMP_pixels(pxl, img);

        /* close file */
        fclose(img);
        return 0;
}

static int_bmp_pixel_t *
readBMP_img_onemalloc(char* filename)
{
        int_bmp_pixel_t  * IntPixelsTab   = NULL;
        FILE *img = NULL;

        /* Openfile */
        img = fopen(filename, "rb");
        assert(img != NULL);

        printf("Read header ...\n");
        readBMP_header(img);
        assert(header != NULL);

        printf("Read pixels ...\n");
        IntPixelsTab = readBMP_pixels_onemalloc(img);

        assert(IntPixelsTab != NULL);

        /* close file */
        fclose(img);

        return IntPixelsTab;
}

static int_bmp_pixel_t **
readBMP_img(char* filename)
{
        int_bmp_pixel_t  ** IntPixelsTab  = NULL;
        FILE *img = NULL;

        /* Openfile */
        img = fopen(filename, "rb");
        assert(img != NULL);

        printf("Read header ...\n");
        readBMP_header(img);
        assert(header != NULL);

        printf("Read pixels ...\n");
        IntPixelsTab = readBMP_pixels(img);

        assert(IntPixelsTab != NULL);

        /* close file */
        fclose(img);

        return IntPixelsTab;
}

/* Accessor header info: Width */
int get_img_width(void)
{
        assert(header != NULL);
        return header->width;
}
int get_img_width_onemalloc(char* filename)
{
        FILE *img = NULL;
        img = fopen(filename, "rb");
        assert(img != NULL);

        printf("Read header ...\n");
        readBMP_header(img);
        assert(header != NULL);
        return header->width;
}

/* Accessor header info: Height */
int get_img_heigh(void)
{
        assert(header != NULL);
        return header->height;
}

/* Accessor header info: Dump infos */
void print_all_header_infos(void)
{
        assert(header != NULL);
        dump_header_bmp(&(header->raw));
}

/* Lecture image */
int_bmp_pixel_t ** Lecture_image( char * filename)
{
        return readBMP_img(filename);
}

/* Lecture image */
int_bmp_pixel_t (*Lecture_image_onemalloc( char * filename))[]
{
        return (int_bmp_pixel_t (*)[header->width])readBMP_img_onemalloc(filename);
}

/* Ecriture image */
int Ecriture_image_onemalloc( int_bmp_pixel_t (*img)[], char* filename)
{
        return writeBMP_img_onemalloc(img, filename);
}
int Ecriture_image( int_bmp_pixel_t ** img, char* filename)
{
        return writeBMP_img(img, filename);
}

void Liberation_image_lue_onemalloc(int_bmp_pixel_t (*img)[])
{
        free(img);
}
void Liberation_image_lue(int_bmp_pixel_t ** img)
{
        int i;
        for(i = 0; i < header->height; i++)
        {
                free(img[i]);
        }
        free(img);
}

void Liberation_finale(void){
        free(header);
}


#ifdef COMPILE_AS_AN_EXECUTABLE
int main( int argc, char * argv[])
{
        int_bmp_pixel_t ** tab;
        /* We need only one argument */
        if( argc != 2)
        {
                fprintf( stderr, "usage: %s <file.bmp>\n", argv[0] );
                return EXIT_FAILURE;
        }

        /* Call read BMP func */
        tab = Lecture_image(argv[1]);
        writeBMP_img(tab, "copie_image");

        return EXIT_SUCCESS;
}
#endif /* COMPILE_AS_AN_EXECUTABLE */
