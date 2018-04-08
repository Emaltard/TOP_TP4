#ifndef __BMP_READER_H__
#define __BMP_READER_H__

struct int_bmp_pixel_s
{
    int Rouge;
    int Bleu;
    int Vert;
};
typedef struct int_bmp_pixel_s int_bmp_pixel_t;

/* Accessor header info: Width */
int get_img_width(void);
int get_img_width_onemalloc(char*);

/* Accessor header info: Height */
int get_img_heigh(void);

/* Accessor header info: Dump infos */
void print_all_header_infos(void);

/* Read BMP img */
int_bmp_pixel_t ** Lecture_image( char *);
int_bmp_pixel_t (*Lecture_image_onemalloc( char *))[];

/* Write BMP img */
int Ecriture_image( int_bmp_pixel_t **, char*);
int Ecriture_image_onemalloc( int_bmp_pixel_t (*img)[], char*);

void Liberation_image_lue(int_bmp_pixel_t ** img);
void Liberation_image_lue_onemalloc(int_bmp_pixel_t (*img)[]);

void Liberation_finale(void);
#endif /* __BMP_READER_H__ */
