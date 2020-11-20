#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include "projet.h"
//Bobillier Laura, Akel Ward


/**
 * @brief  structure d'un point/vecteur
 * 
 */
struct vec {
double x;
double y;
};

/**
 * @brief structure correspondant à une liste de vecteur 
 * 
 */
struct vecset {
struct vec *data;
size_t size;
size_t capacity;
};

/****************Partie 1 : Structures pour la géométrie***********************/

/**
 * @brief calcule le produit scalaire de deux vecteurs
 * 
 * @param v1  le premier vecteur    
 * @param v2 le second vecteur
 * @return double le produit scalaire 
 */
double dot(const struct vec *v1, const struct vec *v2){
  return (v1->x*v2->x)+(v1->y*v2->y);
}

/**
 * @brief fait le produit vectoriel de 2 vecteurs (p1,p2) et (p1,p3)
 * 
 * @param p1 
 * @param p2 
 * @param p3 
 * @return double >0 si angle dans le sens trigonométrique, 0 si collinéaire ,
 *  <0 si dans le sens horaire
 */
double cross(const struct vec *p1,const struct vec *p2, const struct vec *p3){
  double a = (p2->x - p1->x)*(p3->y-p1->y);
  double b = (p2->y-p1->y)*(p3->x-p1->x);
  return a-b;
}



/***********************Partie 2 : Ensemble de points*************************/


/**
 * @brief calcul si les 3 vecteurs forment un tournant à gauche 
 * 
 * @param p1 
 * @param p2 
 * @param p3 
 * @return true si p1,p2,p3 forme un tourne à gauche
 * @return false sinon
 */
bool is_left_turn(const struct vec *p1,const struct vec *p2, const struct vec *p3){
  if (cross(p1,p2,p3)<0){
    return true;
  }
  return false;
}



/**
 * @brief initialise un tableau vide de vecteur
 * 
 * @param self pointeur de la liste
 */
void vecset_create(struct vecset *self){
  struct vec *data = calloc (6, sizeof (const struct vec));
  self ->data = data;
  self->size =0;
  self->capacity = 6;
}

/**
 * @brief détruit un tableau de vecteur
 * 
 * @param self pointeur de structure
 */
void vecset_destroy(struct vecset *self){
  free (self->data);
  self->capacity =0;
  self->size =0;
    
}

/**
 * @brief permet de multiplier la taille d'un tableau dynamique par 2 
 * 
 * @param self pointeur de structure 
 */
void grow_data (struct vecset *self){
  struct vec *data = calloc(self->capacity*2, sizeof(struct vec));  
  memcpy(data,self->data,self->size*sizeof(struct vec));
  free(self->data);
  self->data = data;
  self->capacity= self->capacity*2;
}


/**
 * @brief ajoute un vecteur à la fin d'un tableau de vecteur 
 * 
 * @param self pointeur de l'ensemble de vecteurs
 * @param p le vecteur à ajouter
 * 
 */
void vecset_add(struct vecset *self, struct vec p){

  if (self-> size >= self->capacity){
    grow_data (self);
  }
  self->data [self->size] = p;
  self->size ++;

}

/**
 * @brief création d'un nouveau type de pointeur, le pointeur vers une contions de comparaison 
 * 
 */
typedef int (*comp_func_t)(const struct vec *p1, const struct vec *p2, const void *ctx);

//
/**
 * @brief compare 2 vecteurs en fonction de leurs coordonnées
 * 
 * @param p1 
 * @param p2 
 * @param ctx 
 * @return -1 si pi<p2 ,0 si p1 =p2 et 1 si p1>p2
 */
int comp(const struct vec *p1, const struct vec *p2, const void *ctx){
  if(p1->y<p2->y){
    return -1;
  }
  else if(p1->y==p2->y) {
    if(p1->x==p2->x){
      return 0;
    }
    if(p1->x<p2->x){
      return -1;
    }
  }
  return 1;
}

/**
 * @brief Renvoie le vecteur maximum d'un tableau de vecteurs
 * 
 * @param self 
 * @param func 
 * @param ctx 
 * @return const struct vec* pointeur vers le vecteur max  
 */
const struct vec *vecset_max(const struct vecset *self,comp_func_t func, const void *ctx){
  struct vec *max = &self->data[0];
  for(size_t i = 0 ; i < self->size ; i++){
    if(func (max, &self->data[i], &ctx) == -1 ){
      max = &self->data[i];
    }
  }
  return max;
}

/**
 * @brief Renvoie le vecteur minimum d'un tableau de vecteurs
 * 
 * @param self 
 * @param func 
 * @param ctx 
 * @return const struct vec* pointeur vers le vecteur min
 */
const struct vec *vecset_min(const struct vecset *self,comp_func_t func, const void *ctx){
  struct vec *min = &self->data[0];
  for(size_t i = 0 ; i < self->size ; i++){
    if(func (min, &self->data[i], &ctx) == 1){
      min = &self->data[i];
    }
  }
  return min;
}



/**
 * @brief Échange deux vecteurs dans un tableau de vecteur aux indices k et l 
 * 
 * @param self 
 * @param k 
 * @param l 
 */
static void array_swap (struct vecset  *self, size_t k, size_t l){
  struct vec temp = self->data[k];
  self->data[k] = self->data[l] ;
  self->data[l] =temp;
}

/**
 * @brief Trier un ensemble de vecteurs, tri à bulle ici (pas les bières X)  )
 * 
 * @param self 
 * @param func 
 * @param ctx 
 */
void vecset_sort(struct vecset *self, comp_func_t func, const void *ctx){
  for(size_t i=0; i<self->size-1; ++i){
    for(size_t j =self->size-1; j>i;--j ){
      if(func(&self->data[j],&self->data[j-1],ctx)==-1){
        array_swap(self,j,j-1);
      }
    }
  }
}

/**
 * @brief Insérer élément au début de la pile
 * 
 * @param self 
 * @param p élément à ajouter
 */
void vecset_push(struct vecset *self, struct vec p){
  assert(self);
  vecset_add(self, p);
}

/**
 * @brief supprimer élément au début de la pile
 * 
 * @param self 
 */
void vecset_pop(struct vecset *self){
  assert(self);
  self->size--;
}

/**
 * @brief Renvoie le premier élément de la pile
 * 
 * @param self 
 * @return const struct vec* le premier élément de la pile
 */
const struct vec *vecset_top(const struct vecset *self){
  assert(self);
  return &self->data[self->size-1];
}


/**
 * @brief Renvoie le second élément de la pile
 * 
 * @param self 
 * @return const struct vec* le second élément
 */
const struct vec *vecset_second(const struct vecset *self){
  assert(self);
  return &self->data[self->size-2];
}


/**
 * @brief test si 2 vecteurs sont identiques 
 * 
 * @param v1 
 * @param v2 
 * @return true si les valeurs sont identiques
 * @return false sinon 
 */ 
bool is_equal (const struct vec *v1, const struct vec *v2){
  return (v1->x ==v2->x && v1->y == v2->y);
}

/**************************partie 3  Marche de Jarvis*****************************/
/**
 * @brief calcule une enveloppe convexe selon la méthode de calcul de la marche de jarvis 
 * 
 * @param in ensemble de vecteurs
 * @param out vecteurs composant l'enveloppe convexe
 */
void jarvis_march(const struct vecset *in, struct vecset *out){
  //si le nombre de points est <3 alors il n'est pas possible de calculer une enveloppe convexe 
  if(in->size<3){
    return;
  }
  // if(in->size < 4){
  //   for(size_t i = 0 ; i < in->size-1 ; ++i){
  //     vecset_add(out, in->data[i]);
  //   }
  //   return;
  // }
  int l = 0; 
  for (size_t i = 1; i < in->size; i++) {
    if (in->data[i].x < in->data[l].x) 
      l = i; 
  }
  size_t p = l, q;
  do{
    vecset_add(out, in->data[p]);
    q = (p+1)%in->size;
    for (int i = 0; i < in->size; i++) {      
      if (is_left_turn(&in->data[p], &in->data[i], &in->data[q])){
        q = i; 
      }
    }
    p=q;
  }while(p!=l);
  
}


/*****************************Partie 4 : Parcours de Graham****************************/

/**
 * @brief calcule le carré de la distance qui sépare 2 points 
 * 
 * @param p1 
 * @param p2 
 * @return double le carré de la distance p1-p2
 */
double distance (const struct vec *p1, const struct vec *p2) 
{ 
  return (p1->x - p2->x)*(p1->x - p2->x) + (p1->y - p2->y)*(p1->y - p2->y); 
} 

/**
 * @brief compare les angles formés par p1,p2,p3 
 * 
 * @param p1 
 * @param p2 
 * @param vp0 
 * @return int 
 */
int compare_angle( const struct vec *p1, const struct vec *p2,const void *vp0) 
{ 
  const struct vec *p0 = (struct vec *)vp0; 
   
   // Find orientation 
  int o = cross(p0, p1, p2); 
  if (o == 0) 
    return (distance(p0, p2) >= distance(p0, p1))? 1 : -1; 

  return (o >0 )? -1: 1; 
} 



/**
 * @brief calcule une enveloppe convexe selon la méthode de calcul du parcours de Graham
 * 
 * @param in ensemble de vecteurs
 * @param out vecteurs qui composent l'enveloppe convexe
 */
void graham_scan(const struct vecset * in, struct vecset * out)
{
  assert(in);
  assert(out);
  if(in->size<3){
    return;
  }
  const struct vec *bas = vecset_min(in, &comp, NULL);
  struct vecset *copy = malloc(sizeof(struct vecset));
  vecset_create(copy);

  for (size_t i = 0; i < in->size; ++i) {
    if (comp(&in->data[i], bas, NULL)) {
      vecset_push(copy, in->data[i]);
    }
  }

  vecset_sort(copy, &compare_angle, bas);
  struct vec *first = copy->data;
  vecset_push(out, *bas);
  vecset_push(out, *first);

  const struct vec *top;
  const struct vec *second;

  for (size_t i = 1; i < copy->size; ++i) {
    top = vecset_top(out);
    second = vecset_second(out);
    const struct vec sup = copy->data[i];
    while (out->size >= 2 && is_left_turn(second, top, &sup)) {
      vecset_pop(out);
      top = vecset_top(out);
      second = vecset_second(out);
    }
    vecset_push(out, copy->data[i]);
  }
  vecset_destroy(copy);
  free(copy);
  copy = NULL;
}


/*************************Partie 5 : Enveloppe rapide*****************************/

/**
 * @brief 
 * 
 * @param S 
 * @param X 
 * @param Y 
 * @return struct vecset* 
 */
struct vecset * findhull(struct vecset * S, const struct vec * X, const struct vec * Y)
{
  if (S->size == 0) {
    struct vecset *tab = malloc(sizeof(struct vecset));
    vecset_create(tab);
    return tab;
  }

    // m est le coefficient principal de (X,Y)
    double m = (X->y - Y->y) / (X->x - Y->x);
    // p est l'ordre à l'origine de (X,Y)
    double p = X->y - m * X->x;

    // Calcule le point le plus éloigné de l'ensemble à partir de (X,Y)
    struct vec *P;
    double res = 0;
    for (int i = 0; i < S->size; ++i) {
      if ((fabs(m * S->data[i].x - S->data[i].y + p) / sqrt(pow(m, 2) + 1)) > res) {
        P = &(S->data[i]);
        res = (fabs(m * S->data[i].x - S->data[i].y + p) / sqrt(pow(m, 2) + 1));
      }
    }

    struct vecset *S1 = malloc(sizeof(struct vecset));
    struct vecset *S2 = malloc(sizeof(struct vecset));
    vecset_create(S1);
    vecset_create(S2);

    for (size_t i = 0; i < S->size; i++) {
      if (S->data[i].x != P->x && S->data[i].y != P->y) {
        const struct vec T = S->data[i];
        // Si le point est à gauche de (X,P)
        if (is_left_turn(X, P, &T)) {
          vecset_add(S1, S->data[i]);
        // Si le point est à gauche de (P,Y)
        } else if (is_left_turn(P, Y, &T)) {
          vecset_add(S2, S->data[i]);
        }
      }
    }

    struct vecset *R1 = findhull(S1, X, P);
    struct vecset *R2 = findhull(S2, P, Y);

    struct vecset *R = malloc(sizeof(struct vecset));
    vecset_create(R);

    // Union de R1, P and R2
    for (int i = 0; i < R1->size; ++i) {
      vecset_add(R, R1->data[i]);
    }
    vecset_add(R, *P);
    for (int i = 0; i < R2->size; ++i) {
      vecset_add(R, R2->data[i]);
    }

    // Détruire les vecset crées
    vecset_destroy(S1);
    vecset_destroy(S2);
    vecset_destroy(R1);
    vecset_destroy(R2);
    free(S1);
    free(S2);
    free(R1);
    free(R2);
    S1 = NULL;
    S2 = NULL;
    R1 = NULL;
    R2 = NULL;

    return R;
}


/**
 * @brief calcule une enveloppe convexe selon la méthode de calcul de l'enveloppe rapide 
 * 
 * @param in ensemble de vecteurs
 * @param out vecteurs qui composent l'enveloppe convexe
 */
void quickhull(const struct vecset *in, struct vecset *out)
{
  assert(in);
  assert(out);
  if(in->size<3){
    return;
  }

  // Point le plus à gauche de l'entrée
  const struct vec A = *vecset_min(in, comp, NULL);
  // Point le plus à droite de l'entrée
  const struct vec B = *vecset_max(in, comp, NULL);

  struct vecset *S1 = malloc(sizeof(struct vecset));
  struct vecset *S2 = malloc(sizeof(struct vecset));
  vecset_create(S1);
  vecset_create(S2);

  // Copier de l'entrée dans S (sans les points les plus à gauche et à droite)
  for (size_t i = 0; i < in->size; i++) {
    if (in->data[i].x != A.x && in->data[i].x != B.x && in->data[i].y != A.y && in->data[i].y != B.y)
    {
      const struct vec test = in->data[i];
      if (is_left_turn(&A, &B, &test)) {
        vecset_push(S1, in->data[i]);
      } else {
        vecset_push(S2, in->data[i]);
      }
    }
  }

  struct vecset *R1 = findhull(S1, &A, &B);
  struct vecset *R2 = findhull(S2, &B, &A);

  // Union de A, R1, B and R2
  vecset_add(out, A);
  for (int i = 0; i < R1->size; ++i) {
    vecset_add(out, R1->data[i]);
  }
  vecset_add(out, B);
  for (int i = 0; i < R2->size; ++i) {
    vecset_add(out, R2->data[i]);
  }

  // Détruire les vecset créés
  vecset_destroy(S1);
  vecset_destroy(S2);
  vecset_destroy(R1);
  vecset_destroy(R2);
  free(S1);
  free(S2);
  free(R1);
  free(R2);
  S1 = NULL;
  S2 = NULL;
  R1 = NULL;
  R2 = NULL;
}


/*************************Partie 6 : Pilote***************************/
/**
 * @brief dans le main afin d'utiliser les différentes méthodes de calcul 
 * il suffi juste de décommenter celle que vous voulez utiliser et de commenter les autres 
 * 
 */
#define BUFSIZE 1024
int main() {
  setbuf(stdout, NULL); // avoid buffering in the output
  char buffer[BUFSIZE];
  fgets(buffer, BUFSIZE , stdin);
  struct vecset in;
  vecset_create(&in);
  struct vecset out;
  vecset_create(&out);  
  
  size_t count = strtol(buffer, NULL, 10);
  for (size_t i = 0; i < count; ++i) {
    struct vec p;
    fgets(buffer, BUFSIZE, stdin);
    char *endptr = buffer;
    p.x = strtod(endptr, &endptr);
    p.y = strtod(endptr, &endptr);
    vecset_add(&in,p);
    
  }

  //jarvis_march(&in,&out);
  //graham_scan(&in,&out);
  //quickhull(&in,&out);


  printf("%zu\n", (&out)->size);
  for (size_t i = 0; i < (&out)->size; ++i) {
    printf("%f %f\n", (&out)->data[i].x,(&out)->data[i].y );
  }
  vecset_destroy(&out);
  vecset_destroy(&in);
}

//./hull-generator 10 | ./hull-viewer ./main
// gcc -Wall -std=c99 -O2 -g -o new new.c -lm
