/**
 * \file
 * \brief
 * \author Lectem
 */
#ifndef VECTOR3_H
#define VECTOR3_H

template <class T>
class Vector3
{
    public:
        Vector3(){}
        Vector3(T x,T y,T z):x(x),y(y),z(z) {}
        T x;
        T y;
        T z;

        T dot (Vector3<T> const & v ) const
        {
            return x*v.x+y*v.y+z*v.z;
        }

        Vector3<T> operator+ (Vector3<T> const & v) const
        {
            return Vector3<T>(x+v.x,y+v.y,z+v.z);
        }

        Vector3<T> operator* (Vector3<T> const & v) const
        {
            return Vector3<T>(
                           y*v.z-z*v.y,
                           z*v.x-x*v.z,
                           x*v.y-y*v.x
                           );
        }


        Vector3<T> operator* (T const & c) const
        {
            return Vector3<T>(x*c,y*c,z*c);
        }
};

template<class T>
Vector3<T> operator* (T const c, Vector3<T> const &v)
{
    return Vector3<T>(c*v.x,c*v.y,c*v.z);
}

#endif // VECTOR3_H
