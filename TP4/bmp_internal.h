#ifndef __BMP_INTERNAL_H__
#define __BMP_INTERNAL_H__

struct header_bmp_s {
  int TailleFichier;    // taille totale du fichier
  int Reserves;         // 4 octets reserves usage futur
  int PositionImage;    // position de l'image dans le fichier  
                        // dÃ©placement par rapport au dÃ©but du fichier
  int TailleEntete;     // taille de l'entÃªte d'information 
  int Width;            // largeur de l'image en pixel
  int Height;           // hauteur de l'image en pixel
  short NombrePlan;       // nombre de plans utilises par la carte 
  short NombreBit;        // nombre de bits utilises pour coder un pixel
  int TypeCompres;      // type de compression utilisee (0 => RIEN)
  int TailleImage;      // taille de l'image (hauteur * largeur)
  int ResolutionH;      // resolution horizontale en pixel par metre
  int ResolutionV;      // resolution verticale en pixel par metre
  int NombreCouleur;    // nombre de couleurs => C dans le tableau
  int CouleurImport;    // nombre de couleurs importantes
};
typedef struct header_bmp_s header_bmp_t;

struct bmp_header_infos_s
{
    int width;
    int height;
    header_bmp_t raw;
};
typedef struct bmp_header_infos_s bmp_header_infos_t;

struct char_bmp_pixel_s
{
    char Rouge;
    char Bleu;
    char Vert;
};
typedef struct char_bmp_pixel_s char_bmp_pixel_t;

#endif /* __BMP_INTERNAL_H__ */
