/*8:*/
#line 28 "./common.w"

__device__ double cuda_convert_radian_to_degree(double angle)
{
angle*= RADIAN_TO_DEGREE;
if(angle<0.0)angle+= 360.0;
return angle;
}

/*:8*//*14:*/
#line 46 "./vector.w"

__device__ void cuda_vector_zero(Vector v)
{
v[0]= v[1]= v[2]= 0.0;
v[3]= 1.0;
}

/*:14*//*16:*/
#line 67 "./vector.w"

__device__ void cuda_vector_homogenise(Vector v)
{
if(1.0==v[3])return;
if(0.0!=v[3]){
v[0]/= v[3];
v[1]/= v[3];
v[2]/= v[3];
v[3]= 1.0;
}
}


/*:16*//*18:*/
#line 93 "./vector.w"

__device__ void cuda_vector_copy(Vector u,const Vector v)
{
u[0]= v[0];
u[1]= v[1];
u[2]= v[2];
u[3]= v[3];
}


/*:18*//*20:*/
#line 113 "./vector.w"

__device__ double cuda_vector_magnitude(const Vector v)
{
return sqrt(v[0]*v[0]+
v[1]*v[1]+
v[2]*v[2]);
}


/*:20*//*22:*/
#line 144 "./vector.w"

__device__ void cuda_vector_normalize(const Vector v,Vector r)
{
double m= cuda_vector_magnitude(v);

if(m> 0.0){
r[0]= v[0]/m;
r[1]= v[1]/m;
r[2]= v[2]/m;
}else{
r[0]= 0.0;
r[1]= 0.0;
r[2]= 0.0;
}
}

/*:22*//*24:*/
#line 174 "./vector.w"

__device__ void cuda_vector_difference(const Vector u,const Vector v,Vector r)
{
r[0]= u[0]-v[0];
r[1]= u[1]-v[1];
r[2]= u[2]-v[2];
r[3]= 1.0;
}


/*:24*//*26:*/
#line 195 "./vector.w"

__device__ double cuda_vector_dot(const Vector u,const Vector v)
{
return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

/*:26*//*28:*/
#line 222 "./vector.w"

__device__ void cuda_vector_cross(const Vector u,const Vector v,Vector r)
{
r[0]= (u[1]*v[2]-u[2]*v[1]);
r[1]= (u[2]*v[0]-u[0]*v[2]);
r[2]= (u[0]*v[1]-u[1]*v[0]);
}


/*:28*//*30:*/
#line 250 "./vector.w"

__device__ double cuda_vector_angle_radian(const Vector u,const Vector v)
{
Vector a,b,c= ZERO_VECTOR;
double angle;
cuda_vector_normalize(u,a);
cuda_vector_normalize(v,b);
angle= acos(vector_dot(a,b));
cuda_vector_cross(u,v,c);
if(c[3]<0.0)return(TWICE_PI-angle);
return angle;
}


/*:30*//*32:*/
#line 274 "./vector.w"

__device__ double cuda_vector_angle_degree(const Vector u,const Vector v)
{
return RADIAN_TO_DEGREE*vector_angle_radian(u,v);
}


/*:32*//*34:*/
#line 294 "./vector.w"

__device__ double cuda_vector_distance(const Vector u,const Vector v)
{
double x,y,z;
x= u[0]-v[0];
y= u[1]-v[1];
z= u[2]-v[2];
return sqrt(x*x+y*y+z*z);
}
#line 1 "./matrix.w"


/*:34*//*219:*/
#line 2118 "./csg.w"

__device__ Containment cuda_is_inside_block(const Vector v,const Primitive*p)
{
if(v[0]<p->b.x0||v[0]> p->b.x1||v[1]<p->b.y0||v[1]
> p->b.y1||v[2]<p->b.z0||v[2]> p->b.z1)
return OUTSIDE;
if(v[0]==p->b.x0||v[0]==p->b.x1||v[1]==p->b.y0||v[1]
==p->b.y1||v[2]==p->b.z0||v[2]==p->b.z1)
return SURFACE;
return INSIDE;
}

/*:219*//*221:*/
#line 2161 "./csg.w"

__device__ Containment cuda_is_inside_sphere(const Vector v,const Primitive*p)
{
double delta= sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
if(delta> p->s.radius)return OUTSIDE;
if(delta==p->s.radius)return SURFACE;
return INSIDE;
}


/*:221*//*225:*/
#line 2221 "./csg.w"

__device__ Containment cuda_is_inside_cylinder(const Vector v,const Primitive*p)
{
double delta;
if(v[1]<p->c.y0||v[1]> p->c.y1)return OUTSIDE;
/*226:*/
#line 2238 "./csg.w"

delta= sqrt(v[0]*v[0]+v[2]*v[2]);

/*:226*/
#line 2226 "./csg.w"
;
if(delta> p->c.radius)return OUTSIDE;
if(v[1]==p->c.y0||v[1]==p->c.y1||delta==
p->c.radius)
return SURFACE;
return INSIDE;
}

/*:225*//*239:*/
#line 2376 "./csg.w"

__device__ Containment cuda_is_inside_torus(const Vector v,const Primitive*p)
{
double gamma,gamma_deg,tau,tau_deg,delta,radial;
Vector tube_center,from_tube_center_to_v,temp;
/*232:*/
#line 2331 "./csg.w"

delta= sqrt(v[0]*v[0]+v[2]*v[2]);

/*:232*/
#line 2381 "./csg.w"
;
if(delta<p->t.r0||delta> p->t.r1)return OUTSIDE;

/*240:*/
#line 2392 "./csg.w"

cuda_vector_copy(temp,v);
temp[1]= 0.0;
gamma= cuda_vector_angle_radian(positive_xaxis_unit_vector,temp);
gamma_deg= cuda_convert_radian_to_degree(gamma);

/*:240*/
#line 2384 "./csg.w"
;
if(angle_outside_range(gamma_deg,p->t.phi_start,p->t.phi_end))
return OUTSIDE;
/*241:*/
#line 2398 "./csg.w"

/*235:*/
#line 2356 "./csg.w"

tube_center[0]= p->t.major*cos(gamma);
tube_center[1]= 0.0;
tube_center[2]= p->t.major*sin(gamma);

/*:235*/
#line 2399 "./csg.w"
;
/*242:*/
#line 2405 "./csg.w"

cuda_vector_difference(v,tube_center,from_tube_center_to_v);
radial= cuda_vector_magnitude(from_tube_center_to_v);

/*:242*/
#line 2400 "./csg.w"
;
if(radial> p->t.minor)return OUTSIDE;
/*243:*/
#line 2409 "./csg.w"

tau= cuda_vector_angle_radian(tube_center,from_tube_center_to_v);
tau_deg= cuda_convert_radian_to_degree(tau);

/*:243*/
#line 2402 "./csg.w"
;
if(angle_outside_range(tau,p->t.theta_start,p->t.theta_end))return OUTSIDE;

/*:241*/
#line 2387 "./csg.w"
;
/*238:*/
#line 2369 "./csg.w"

if(radial==p->t.minor)return SURFACE;
if(p->t.phi<360.0&&(gamma==p->t.phi_start||gamma==p->t.phi_end))
return SURFACE;
if(p->t.theta<360.0&&(tau==p->t.theta_start||tau==p->t.theta_end))
return SURFACE;

/*:238*/
#line 2388 "./csg.w"
;
return INSIDE;
}

/*:239*//*317:*/
#line 41 "./bstack.w"

#define MAX_BOOLEAN_STACK_SIZE 1024
typedef struct{
int tos,size;
bool v[MAX_BOOLEAN_STACK_SIZE];
}boolean_stack;
__device__ bool boolean_stack_init(boolean_stack*s)
{
if(NULL==s)return false;
s->tos= 0;
s->size= MAX_BOOLEAN_STACK_SIZE;
return true;
}
__device__ bool boolean_stack_push(boolean_stack*s,bool v)
{
if(s->tos==s->size)return false;
s->v[s->tos++]= v;
return true;
}
__device__ bool boolean_stack_pop(boolean_stack*s,bool*v)
{
if(0==s->tos)return false;
*v= s->v[--s->tos];
return true;
}

#line 1 "./postfix.w"


/*:317*//*320:*/
#line 29 "./postfix.w"

__device__ bool cuda_is_inside_primitive(Vector v,Primitive*p)
{
Containment c;
switch(p->type){
case BLOCK:c= cuda_is_inside_block(v,p);break;
case SPHERE:c= cuda_is_inside_sphere(v,p);break;
case CYLINDER:c= cuda_is_inside_cylinder(v,p);break;
case TORUS:c= cuda_is_inside_torus(v,p);break;
default:c= INVALID;
}
if(INSIDE==c||SURFACE==c)return true;
return false;
}

/*:320*/
