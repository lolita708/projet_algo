#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "projet.h"

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
  self ->data = calloc (6, sizeof (const struct vec));
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

static void array_swap (struct vecset  *self, size_t k, size_t l){
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
    printf("(%f , %f) \n", self->data[i].x, self->data[i].y);
  }
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
  // const struct vec *start = vecset_min (in, comp, NULL);
  // vecset_add(out, *start);
  // struct vec prevVec = *start;
  // struct vec *next;
  
  // do{
  //   for(size_t i =0; i< in->size; ++i){
      
  //   }
    
  // }while(!is_equal(start,&prevVec));
}


////Partie 4 : Parcours de Graham/////

void graham_scan(const struct vecset *in, struct vecset *out){
  //Le point le plus bas
  const struct vec *bas = vecset_min (in, comp, NULL);
  // Copier l'entrée dans copy
  struct vecset *copy = malloc(sizeof(struct vecset));
  vecset_create(copy);
  
  for (size_t i = 0; i < in->size; ++i) {
    if(comp(&in->data[i], bas, NULL)){
      vecset_push(copy, in->data[i]);
    }
  }
  


}


//////Partie 5 : Enveloppe rapide/////

void quickhull(const struct vecset *in, struct vecset *out);


/////Partie 6 : Pilote//////
#define BUFSIZE 1024
int main() {
  struct vecset *in;
  struct vecset *out;
  vecset_create(in);
  vecset_create(out);
  setbuf(stdout, NULL); // avoid buffering in the output
  char buffer[BUFSIZE];
  fgets(buffer, BUFSIZE , stdin);
  size_t count = strtol(buffer, NULL, 10);
  for (size_t i = 0; i < count; ++i) {
    struct vec p;
    fgets(buffer, BUFSIZE, stdin);
    char *endptr = buffer;
    p.x = strtod(endptr, &endptr);
    p.y = strtod(endptr, &endptr);
    vecset_add(in,p);
    // then do something with p
  }

  jarvis_march(in, out);
  

  return 0;
}

//./hull-generator 10 | ./hull-viewer ./main