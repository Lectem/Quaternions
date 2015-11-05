#ifndef QUATERNIONS_H_
#define QUATERNIONS_H_

#include <cmath>
#include <cfloat> //Nécessaire pour Epsilon
#include "Vector3.h"



template <class T=float>
class Quaternion
{
public:
    T w,x,y,z;
    Quaternion() {}
    Quaternion(T w,T x,T y,T z):w(w),x(x),y(y),z(z) {}
    Quaternion(T w,const Vector3<T> &v):w(w),x(v.x),y(v.y),z(v.z) {}

    T dot(Quaternion<T> const & q) const;
    T norm2() const;
    T norm() const;
    Quaternion<T> getConjugate() const
    {
        return Quaternion<T>(w,-x,-y,-z);
    }

    Quaternion<T> operator+ (const Quaternion<T> &rq) const;
    Quaternion<T>& operator+= (const Quaternion<T> &rq);

    Quaternion<T> operator- (const Quaternion<T> &rq) const;
    Quaternion<T>& operator-= (const Quaternion<T> &rq);

    Quaternion<T> operator* (const T c) const;
    Quaternion<T> operator* (const Quaternion<T> &rq) const;
    Quaternion<T>& operator*= (const Quaternion<T> &rq);

    Vector3<T> operator* (const Vector3<T> &vec) const;

    Quaternion<T> operator/ (const T c) const;

    void fromVector3(const Vector3<T> &vec);
    void FromAxisAngle(const Vector3<T> &v,T theta);
    void FromEuler(T roll,T pitch,T yaw);

    void toRotMat4x4(T mat[16]) const;
    void toRotMat3x3(T mat[9]) const;

    bool is_unit_length() const;
    void normalize();

    static Quaternion<T> slerp(Quaternion<T> const & q1,Quaternion<T> const & q2,T t);
    static Quaternion<T> nlerp(Quaternion<T> const & q1,Quaternion<T> const & q2,T t);
    static Quaternion<T> squad( Quaternion<T> const & q1,Quaternion<T> const & q2,
                                Quaternion<T> const & a,Quaternion<T> const & b,T t);

    void setEpsilon(T v)
    {
        Epsilon=v;
    }
    T getEpsilon()
    {
        return Epsilon;
    }

private:
    void normalize_with_check();
    static T Epsilon;
};

template <> float Quaternion<float>::Epsilon=10*FLT_EPSILON;
template <> double Quaternion<double>::Epsilon=10*DBL_EPSILON;

/**< Renvoie le produit scalaire avec q */
template <class T>
T Quaternion<T>::dot(Quaternion<T> const & q) const
{
    return w * q.w + x * q.x + y * q.y + z * q.z;
}

/**< Renvoie la norme au carré */
template <class T>
T Quaternion<T>::norm2() const
{
    return w * w + x * x + y * y + z * z;
}

/**< Renvoie la norme */
template <class T>
T Quaternion<T>::norm() const
{
    return sqrt(norm2());
}

/**< Le quaternion est-il unitaire ? */
template <class T>
bool Quaternion<T>::is_unit_length() const
{
    return fabs(norm2() - 1.0f) <= Epsilon;
}

/**< Normalise un quaternion. Si il est presque normalisé, on ne fait pas l'opération pour gagner du temps */
template <class T>
void Quaternion<T>::normalize_with_check()
{
    T mag2 = w * w + x * x + y * y + z * z;
    /*On vérifie si le quaternion est presque normalisé, et si il n'est pas proche de 0*/
    if (!is_unit_length())
    {

        T mag = sqrt(mag2);
        w /= mag;
        x /= mag;
        y /= mag;
        z /= mag;
    }
}

/**< Normalise un quaternion. Assez lent à cause de sqrt et des 4 divisions.*/
template <class T>
void Quaternion<T>::normalize()
{
    T mag = sqrt(w * w + x * x + y * y + z * z);
    w /= mag;
    x /= mag;
    y /= mag;
    z /= mag;
}


/**< Addition */
template <class T>
Quaternion<T> Quaternion<T>::operator+ (const Quaternion<T> &rq) const
{
    return Quaternion<T>(w+rq.w,
                         x+rq.x,
                         y+rq.y,
                         z+rq.z);
}
template <class T>
Quaternion<T>& Quaternion<T>::operator+= (const Quaternion<T> &rq)
{
    w+=rq.w;
    x+=rq.x;
    y+=rq.y;
    z+=rq.z;
    return *this;
}

/**< Soustraction */
template <class T>
Quaternion<T> Quaternion<T>::operator- (const Quaternion<T> &rq) const
{
    return Quaternion<T>(w-rq.w,
                         x-rq.x,
                         y-rq.y,
                         z-rq.z);
}
template <class T>
Quaternion<T>& Quaternion<T>::operator-= (const Quaternion<T> &rq)
{
    w-=rq.w;
    x-=rq.x;
    y-=rq.y;
    z-=rq.z;
    return *this;
}

/**< Multiplication. Attention celle-ci n'est pas commutative*/
template <class T>
Quaternion<T> Quaternion<T>::operator* (const Quaternion<T> &rq) const
{
    return Quaternion<T>(
               w * rq.w - x * rq.x - y * rq.y - z * rq.z,
               w * rq.x + x * rq.w + y * rq.z - z * rq.y,
               w * rq.y + y * rq.w + z * rq.x - x * rq.z,
               w * rq.z + z * rq.w + x * rq.y - y * rq.x);
}
template <class T>
Quaternion<T> & Quaternion<T>::operator *= (Quaternion<T> const & rq)
{
    T wt = w * rq.w - x * rq.x - y * rq.y - z * rq.z;
    T xt = w * rq.x + x * rq.w + y * rq.z - z * rq.y;
    T yt = w * rq.y + y * rq.w + z * rq.x - x * rq.z;
    T zt = w * rq.z + z * rq.w + x * rq.y - y * rq.x;

    w = wt;
    x = xt;
    y = yt;
    z = zt;

    return(*this);
}

/**< Multiplication avec un scalaire*/
template <class T>
Quaternion<T> Quaternion<T>::operator* (T const c) const
{
    return Quaternion<T>(w*c,x*c,y*c,z*c);
}

/**< Division par un scalaire*/
template <class T>
Quaternion<T> Quaternion<T>::operator/ (T const c) const
{
    return Quaternion<T>(w/c,x/c,y/c,z/c);
}


template <class T>
void Quaternion<T>::fromVector3(const Vector3<T> &vec)
{
    x = vec.x;
    y = vec.y;
    z = vec.z;
    w =0.0;
}


/**< rotation gloutonne d'un vecteur 3D par un quaternion. /!\ N'utiliser que des quaternions unitaires.*/
//template <class T>
//Vector3<T> Quaternion<T>::operator* (const Vector3<T> &vec) const
//{
//    //normalize_with_check();
//    Quaternion<T> resQuat(0,vec) ;
//    resQuat *= getConjugate();
//    resQuat = (*this) * resQuat;
//    return (Vector3<T>(resQuat.x, resQuat.y, resQuat.z));
//}

/**< rotation d'un vecteur 3D par un quaternion. /!\ N'utiliser que des quaternions unitaires.*/
template <class T>
Vector3<T> Quaternion<T>::operator* (const Vector3<T> &v) const
{
    //normalize_with_check();
    Vector3<T> r(x,y,z);
    return v+(r+r)*(r*v+w*v);
}

/**< Crée le quaternion correspondant à la rotation d'axe dirigé par (x,y,z) et d'angle theta.
     v doit être unitaire */
template <class T>
void Quaternion<T>::FromAxisAngle(const Vector3<T> &v,T theta)
{
    theta*=0.5;
    T sin_theta=sin(theta);
    w=cos(theta);
    x=v.x*sin_theta;
    y=v.y*sin_theta;
    z=v.z*sin_theta;
    normalize_with_check();
}


/**<    Crée le quaternion correspondant à la rotation d'euler
        Les rotations sont faites dans l'ordre yaw (axe z) > pitch (axe y) > roll (axe x) */
template <class T>
void Quaternion<T>::FromEuler( T roll,T pitch,T yaw)
{
    /*  Version rapide, en réalité est équivalent à créer les 3 quaternions de rotation
        et de les multiplier ensemble                                                   */

    T y = yaw   / 2.0;
    T p = pitch / 2.0;
    T r = roll  / 2.0;

    T siny = sin(y);
    T sinp = sin(p);
    T sinr = sin(r);
    T cosy = cos(y);
    T cosp = cos(p);
    T cosr = cos(r);

    this->w = cosr * cosp * cosy + sinr * sinp * siny;
    this->x = sinr * cosp * cosy - cosr * sinp * siny;
    this->y = cosr * sinp * cosy + sinr * cosp * siny;
    this->z = cosr * cosp * siny - sinr * sinp * cosy;
}


/**< La formule utilisée ne fonctionne que pour des quaternions unitaires.
La matrice de retour s'utilise de la même façon qu'une matrice de rotation standard sous OpenGL.
NB : Attention, en mémoire les vecteurs sont continus, ce qui correspond à des vecteurs à l'horizontale
en notation mathématique */
template <class T>
void Quaternion<T>::toRotMat4x4(T mat[16]) const
{
    T xx      = x * x;
    T xy      = x * y;
    T xz      = x * z;
    T xw      = x * w;

    T yy      = y * y;
    T yz      = y * z;
    T yw      = y * w;

    T zz      = z * z;
    T zw      = z * w;

    mat[0]  = 1 - 2 * ( yy + zz );
    mat[1]  =     2 * ( xy + zw );
    mat[2]  =     2 * ( xz - yw );


    mat[4]  =     2 * ( xy - zw );
    mat[5]  = 1 - 2 * ( xx + zz );
    mat[6]  =     2 * ( yz + xw );

    mat[8]  =     2 * ( xz + yw );
    mat[9]  =     2 * ( yz - xw );
    mat[10] = 1 - 2 * ( xx + yy );

    mat[3]  = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0;
    mat[15] = 1;
}


/**< La formule utilisée ne fonctionne que pour des quaternions unitaires.
La matrice de retour s'utilise de la même façon qu'une matrice de rotation standard sous OpenGL.
NB : Attention, en mémoire les vecteurs sont contigus, ce qui correspond à des vecteurs à l'horizontale
en notation mathématique */
template <class T>
void Quaternion<T>::toRotMat3x3(T mat[9]) const
{
    T xx      = x * x;
    T xy      = x * y;
    T xz      = x * z;
    T xw      = x * w;

    T yy      = y * y;
    T yz      = y * z;
    T yw      = y * w;

    T zz      = z * z;
    T zw      = z * w;

    mat[0]  = 1 - 2 * ( yy + zz );
    mat[1]  =     2 * ( xy + zw );
    mat[2]  =     2 * ( xz - yw );


    mat[3]  =     2 * ( xy - zw );
    mat[4]  = 1 - 2 * ( xx + zz );
    mat[5]  =     2 * ( yz + xw );

    mat[6]  =     2 * ( xz + yw );
    mat[7]  =     2 * ( yz - xw );
    mat[8] = 1 - 2 * ( xx + yy );
}


template <class T>
Quaternion<T> Quaternion<T>::slerp(Quaternion<T> const & q1,Quaternion<T> const & q2,T t)
{
    T theta = acos(q1.dot(q2)),mult1,mult2,sintheta;
    // lorsque theta est très petit,on fait une LERP
    // cela permet d'éviter la division par 0
    if (theta > 0.000001)
    {
        mult1 = sin( (1-t)*theta );;
        mult2 = sin( t*theta );
        sintheta=sin(theta);
    }
    else
    {
        mult1 = 1 - t;
        mult2 = t;
        sintheta=1;
    }
    return (q1*mult1 + q2*mult2)/sintheta;
}

/* Effectue une NLERP entre deux quaternions unitaires */
template <class T>
Quaternion<T> Quaternion<T>::nlerp(Quaternion<T> const & q1,Quaternion<T> const & q2,T t)
{
    Quaternion q(q1*(1-t) + q2*t);
    q.normalize();
    return q;
}


/*  Effectue une interpolation SQUAD entre deux quaternions unitaires
    a et b sont deux quaternions intermédiaires qui dirigent l'interpolation*/
template <class T>
Quaternion<T> Quaternion<T>::
    squad(Quaternion<T> const & q1,Quaternion<T> const & q2,
          Quaternion<T> const & a,Quaternion<T> const & b,T t)
{
    return slerp(slerp(q1,q2,t),slerp(a,b,t),2*t*(1-t));
}


#endif
