#ifndef H_INTERPOLATOR
#define H_INTERPOLATOR

template<class T> class Interpolator {
public:
    /*
        Linear interpolation
        target  - the target point, 0.0 - 1.0
        v       - a pointer to an array of size 2 containg the two values    
    */

    inline static T Linear(float target, T *v) {
        return (T) (target * (v[1]) + (T(1.0f) - target) * (v[0]));
    }

    /*
        BiLinear interpolation, linear interpolation in 2D
        target  - a 2D point (X,Y)
        v       - an array of size 4 containg target values 
                v[3]------v[2]
                 |          |
                 |          |
                 |          |
                 |          |
                v[0]------v[1]

    */
    inline static T Bilinear(float *target, T *v) {
        T v_prime[2] = {
                Linear(target[1], &(v[0])),
                Linear(target[1], &(v[2]))
        };

        return Linear(target[0], v_prime);

    }

    /* 
        TriLinear interpolation, linear interpolation in 2D
        target  - a 3D point (X,Y,Z)
        v       - an array of size 8 containg the target values of the 8 conors
                v[3]------v[2]
              /  |       /  |
             /---|----- /   |
            v[7] |     v[6] |
             |   |     |    |
             |  v[0]------v[1]
             | /       | /
            v[4]-------v[5]      
    */

    inline static T Trilinear(float *target, T *v) {
        T v_prime[2] = {
                Bilinear(&(target[0]), &(v[0])),
                Bilinear(&(target[1]), &(v[4]))
        };

        return Linear(target[2], v_prime);
    }
};

#endif
