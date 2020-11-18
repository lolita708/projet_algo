#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "projet.h"
//test
// /**
//  * @brief calcule le produit scalaire de deux vecteurs
//  * 
//  * @param v1  le premier vecteur    
//  * @param v2 le second vecteur
//  * @return double le produit scalaire 
//  */
double dot(const struct vec *v1, const struct vec *v2){
  return (v1->x*v2->x)+(v1->y*v2->y);
}

/**
 * @brief fait le produit vectoriel de 2 vecteurs (p1,p2) et (p1,p3)
 * 
 * @param p1 
 * @param p2 
 * @param p3 
 * @return double >0 si angle dans le sens trigonomÃ©trique, 0 si collinÃ©raire , <0 si dans l'inverse du sens trigo
 */
double cross(const struct vec *p1,const struct vec *p2, const struct vec *p3){
  double a = (p2->x - p1->x)*(p3->y-p1->y);
  double b = (p2->y-p1->y)*(p3->x-p1->x);
  return a-b;
}

/**
 * @brief 
 * 
 * @param p1 
 * @param p2 
 * @param p3 
 * @return true si Ã§a tourne Ã  gauche
 * @return false sinon
 */
bool is_left_turn(const struct vec *p1,const struct vec *p2, const struct vec *p3){
  if (cross(p1,p2,p3)<0){
    return true;
  }
  return false;
}
///////////vectset
void vecset_create(struct vecset *self){
  struct vec *data = calloc (6, sizeof (const struct vec));
  self ->data = data;
  self->size =0;
  self->capacity = 6;
}

void vecset_destroy(struct vecset *self){
  free (self->data);
  self->capacity =0;
  self->size =0;
    
}

void grow_data (struct vecset *self){
  struct vec *data = calloc(self->capacity*2, sizeof(struct vec));  
  memcpy(data,self->data,self->size*sizeof(struct vec));
  free(self->data);
  self->data = data;
  self->capacity= self->capacity*2;
}



void vecset_add(struct vecset *self, struct vec p){

  if (self-> size >= self->capacity){
    grow_data (self);
  }
  self->data [self->size] = p;
  self->size ++;

}

typedef int (*comp_func_t)(const struct vec *p1, const struct vec *p2, const void *ctx);

//-1 si pi<p2
//0 si =
//1 si p1>p2

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


const struct vec *vecset_max(const struct vecset *self,comp_func_t func, const void *ctx){
  struct vec *max = &self->data[0];
  for(size_t i = 0 ; i < self->size ; i++){
    if(func (max, &self->data[i], &ctx) == -1 ){
      max = &self->data[i];
    }
  }
  return max;
}

const struct vec *vecset_min(const struct vecset *self,comp_func_t func, const void *ctx){
  struct vec *min = &self->data[0];
  for(size_t i = 0 ; i < self->size ; i++){
    if(func (min, &self->data[i], &ctx) == 1){
      min = &self->data[i];
    }
  }
  return min;
}

static void array_swap (struct vecset *self, size_t k, size_t l){
  struct vec temp = self->data[k];
  self->data[k] = self->data[l] ;
  self->data[l] =temp;
}

size_t array_partition(struct vecset *self, comp_func_t func, const void *ctx, size_t i, size_t j) {
  size_t pivot_index = i;
  const struct vec pivot = self->data[pivot_index];
  array_swap(self,pivot_index, j);
  size_t l = i;
  for (size_t k=i;k<j;++k){
    if(func(&self->data[k],&pivot,ctx) == -1){
      array_swap(self,k,l);
      ++l;
    }
  }
  array_swap(self, l, j);
  return l;
}
static void vecset_sort_partial(struct vecset *self, comp_func_t func, const void *ctx,size_t i, size_t j) {
  if(i<j){
    size_t p = array_partition(self,func,ctx, i,j);
    vecset_sort_partial(self,func,ctx, i, p-1);
    vecset_sort_partial(self,func,ctx, p+1, j);    
  }
}

void vecset_sort(struct vecset *self, comp_func_t func, const void *ctx){
  vecset_sort_partial(self, func, ctx, 0, self->size-1);
}

void vecset_push(struct vecset *self, struct vec p){
  assert(self);
  vecset_add(self, p);
}
void vecset_pop(struct vecset *self){
  assert(self);
  self->size--;
}
const struct vec *vecset_top(const struct vecset *self){
  assert(self);
  return &self->data[self->size-1];
}

const struct vec *vecset_second(const struct vecset *self){
  assert(self);
  return &self->data[self->size-2];
}

void afficher_vecset(const struct vecset *self){
  for(size_t i =0; i<self->size;++i){
    printf("%f %f\n", self->data[i].x, self->data[i].y);
  }
}

static int compare_angle(const struct vec *p1, const struct vec *p2, const void *ctx)
{
    assert(ctx);
    assert(p1);
    assert(p2);
    const struct vec * ref = ctx;

    double X1 = p1->x - ref->x;
    double Y2 = p1->y - ref->y;

    double X3 = p2->x - ref->x;
    double Y4 = p2->y - ref->y;

    double a1 = atan2(X1, Y2);
    double a2 = atan2(X3, Y4);

    if (a1 > a2) {
        return 1;
    }
    if (a1 < a2) {
        return -1;
    }
    return comp(p1, p2, ctx);
}

///////////partie 3///// Marche de Jarvis

bool is_equal (const struct vec *v1, const struct vec *v2){
  return (v1->x ==v2->x && v1->y == v2->y);
}

void jarvis_march(const struct vecset *in, struct vecset *out){
  if(in->size < 4){
    for(size_t i = 0 ; i < in->size-1 ; ++i){
      vecset_add(out, in->data[i]);
    }
    return;
  }
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


////Partie 4 : Parcours de Graham/////

void graham_scan(const struct vecset *in, struct vecset *out){
  assert(in);
  assert(out);
  //Le point le plus bas
  struct vec *bas = vecset_min (in, comp, NULL);
  // Copier l'entrée dans copy
  struct vecset copy ;
  vecset_create(&copy);
  
  for (size_t i = 0; i < in->size; ++i) {
    if(comp(&in->data[i], bas, NULL)){
      vecset_push(&copy, in->data[i]);
    }
  }
  //Sort
  vecset_sort(&copy, &compare_angle, bas);

  //Premier élément
  struct vec *first = (&copy)->data;

  //push le point le plus bas et le premier point dans out
  vecset_push(out, *bas);
  vecset_push(out, *first);

  const struct vec *top;
  const struct vec *second;

  for (size_t i = 1; i <(&copy)->size; ++i) {
    top = vecset_top(out);
    second = vecset_second(out);
    const struct vec sup = (&copy)->data[i];
    while (out->size >= 2 && is_left_turn(second, top, &sup)){
      vecset_pop(out);
      top = vecset_top(out);
      second = vecset_second(out);
    }
    vecset_push(out, (&copy)->data[i]);
  }

  //Détruire la copie de in
  vecset_destroy(&copy);
}


//////Partie 5 : Enveloppe rapide/////

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
    for (size_t i = 0; i < S->size; ++i) {
      if ((fabs(m * S->data[i].x - S->data[i].y + p) / sqrt(pow(m, 2) + 1)) > res) {
        P = &S->data[i];
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
    for (size_t i = 0; i < R1->size; ++i) {
        vecset_add(R, R1->data[i]);
    }
    vecset_add(R, *P);
    for (size_t i = 0; i < R2->size; ++i) {
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



void quickhull(const struct vecset *in, struct vecset *out)
{
    assert(in);
    assert(out);

    // Point le plus à gauche de l'entrée
    const struct vec A = *vecset_min(in, compare_x, NULL);
    // Point le plus à droite de l'entrée
    const struct vec B = *vecset_max(in, compare_x, NULL);

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
}


/////Partie 6 : Pilote//////
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

  jarvis_march(&in,&out);
  
  //char buffer2[BUFSIZE];
  printf("%zu\n", (&out)->size);
  for (size_t i = 0; i < (&out)->size; ++i) {
    // fgets(buffer2, BUFSIZE,(&out)->data[i]);

    // char *endptr = buffer;
    // double x = strtod(endptr, &endptr);
    // double y = strtod(endptr, &endptr);

    printf("%f %f\n", (&out)->data[i].x,(&out)->data[i].y );
  }
  vecset_destroy(&out);
}

//./hull-generator 10 | ./hull-viewer ./main
// gcc -Wall -std=c99 -O2 -g -o main main.c
// valgrind --leak-check=full ./main



