struct vec {
double x;
double y;
};

struct vecset {
struct vec *data;
size_t size;
size_t capacity;
};
