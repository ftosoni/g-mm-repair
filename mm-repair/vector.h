/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * Basic definitions for 1-dimensional dense vectors
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

// one dimensional dense uncompressed vector 
typedef struct {
  matval *v;
  size_t size;
} vector;


vector *vector_create();
vector *vector_create_value(int n,matval v);
void vector_set_zero(vector *v,int n);
vector *vector_split(vector *v,int n);
xmatval vector_normalize(vector *v);
void vector_destroy(vector *v);


// ---------- vectors -----------------

// return a pointer to an empty vector object
vector *vector_create()
{
  vector *w = (vector *) malloc(sizeof(vector));
  if(w==NULL) die("malloc failed");
  w->v = NULL;
  w->size=0;
  return w;
}

// return a pointer to a vector of a given size with constant value
vector *vector_create_value(int n, matval v)
{
  vector *w = vector_create();
  assert(w!=NULL);
  w->size=n;
  w->v = (matval *) malloc(n*sizeof(matval));
  if(w->v==NULL) die("malloc failed");
  for(int i=0;i<n;i++) 
    w->v[i]=v;
  return w;
}


// split a vector into an array of n subvectors
// since these subvectors share the same storage
// and do not have independent life, they are not  
// vector pointers, but just elements of an array 
vector *vector_split(vector *v, int n)
{
  assert(n>1);
  assert(v!=NULL);
  vector *w = (vector *) malloc(n*sizeof *w);
  if(w==NULL) die("Out of memory");
  int maxblock = (v->size+n-1)/n;
  assert(maxblock>=1);
  
  int remaining = v->size;
  for(int i=0;i<n;i++) {
    assert(remaining>0);
    int r = (remaining>maxblock? maxblock : remaining);
    assert(r>0);
    // init subvector
    w[i].v = v->v + (v->size-remaining);
    w[i].size = r;
    remaining -= r;
  }
  assert(remaining==0);
  return w;
}



// infinity norm normalization, return the norm before normalization
xmatval vector_normalize(vector *w)
{
  xmatval norm = 0;
  for(int i=0;i<w->size;i++) {
    xmatval t = w->v[i]>0 ? w->v[i] : -w->v[i];
    if(t>norm) norm=t;
  }
  assert(norm>=0);
  if(norm>0) 
    for(int i=0;i<w->size;i++)
      w->v[i] = w->v[i]/norm;
  return norm;
}  

// set v = v + t w
void vector_update(vector *v, xmatval t, vector *w)
{
  assert(v->v!=NULL && w->v!=NULL);
  assert(v->size==w->size);
  for(int i=0;i<v->size;i++)
    v->v[i] += t*w->v[i];
}

// set v = tv
void vector_scalar_update(vector *v, xmatval t)
{
  assert(v->v!=NULL);
  for(int i=0;i<v->size;i++)
    v->v[i] *= t;
}

// norm2 squared computation
xmatval vector_norm2sq(vector *w)
{
  xmatval norm = 0;
  for(int i=0;i<w->size;i++) {
    norm +=   w->v[i] * w->v[i];
  }
  return norm;  
}

// scalar product of two vectors 
xmatval vector_scalar_prod(vector *v, vector *w)
{
  assert(v->size==w->size);
  xmatval sp = 0;
  for(int i=0;i<w->size;i++) {
    sp +=   w->v[i] * w->v[i];
  }
  return sp;  
}
// set an existing vector to 0
void vector_set_zero(vector *v, int dim)
{
  if(v->size != dim) {
    v->v = (matval *) realloc(v->v,dim*sizeof(matval));
    if(v->v==NULL) die("Realloc failed");
    v->size = dim;
  }
  for(int i=0;i<v->size;i++) v->v[i]=0;  
}

void vector_destroy(vector *v)
{
  if(v->v!=NULL) free(v->v);
  free(v);
}


