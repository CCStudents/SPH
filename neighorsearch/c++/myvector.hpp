// User のベクトル定義
// This is slightly modified version of vector header file
#ifndef  STARLAB_VECTOR_H
#define  STARLAB_VECTOR_H

#ifdef PERIODIC
  inline real readjust_r ( real x )
  {
    return fmod(x + 2.5, 1.0) - 0.5;
  }
#else
  inline real readjust_r ( real x )
  {
    return x;
  }
#endif

//引き出すときはvector でもmyvectorでもオーケー
#define vector myvector

//3次元を固定
const int ndim = 3;

class vector{
  private:
    real element[3];
  public:
    vector ( real c = 0 ) // default: 0 で初期化
    {
      element[0] = element[1] = element[2] = c;
    }

    vector ( real x, real y, real z )
    {
      element[0] = x; element[1] = y; element[2] = z;
    }

    // i番目の要素を引き出す
    real & operator [] ( int i ) {return element[i];}

    inline void print ()
    {
      std::cout << element[0] << " "
                << element[1] << " "
			          << element[2] << "\n";
    }
    #ifdef PERIODIC
      vector readjust ()
      {
      	return vector(readjust_r(element[0]), readjust_r(element[1]), readjust_r(element[2]));
      }
    #else
      vector readjust ()
      {
    	   return vector(*this);
      }
    #endif

    // 逆符号
    vector operator - ()
    {
      return vector(-element[0], -element[1], -element[2]);
    }

    // 内積
    real operator * ( const vector& b )
	  {
      return  element[0] * b.element[0]
            + element[1] * b.element[1]
  		      + element[2] * b.element[2];
    }

    // 外積
    vector operator ^ ( const vector &b )
    {
      return vector(element[1] * b.element[2] - element[2] * b.element[1],
  		              element[2] * b.element[0] - element[0] * b.element[2],
  		              element[0] * b.element[1] - element[1] * b.element[0]);
    }

    // 足し算/引き算
    vector operator + ( const vector &b )
	  {
      return vector(element[0] + b.element[0],
                    element[1] + b.element[1],
                    element[2] + b.element[2]);
    }
    vector operator - ( const vector &b )
	  {
      return vector(element[0] - b.element[0],
              	 	  element[1] - b.element[1],
              		  element[2] - b.element[2]);
    }

    friend vector operator + (real, const vector & );
    friend vector operator + (const vector &, real);

    // スカラー倍
    friend vector operator * (real, const vector & );
    friend vector operator * (const vector &, real);
    friend vector operator / (const vector &, real);

    // +=, -=, *=, /= 演算
    vector & operator += ( const vector& b )
	  {
      element[0] += b.element[0];
	    element[1] += b.element[1];
	    element[2] += b.element[2];
	    return *this;
    }
	  vector & operator -= ( const vector& b )
	  {
      element[0] -= b.element[0];
	    element[1] -= b.element[1];
	    element[2] -= b.element[2];
	    return *this;
    }
	  vector & operator *= ( const real b )
	  {
      element[0] *= b;
	    element[1] *= b;
	    element[2] *= b;
	    return *this;
    }
	  vector & operator /= ( const real b )
	  {
      register real binv = 1.0 / b;
      element[0] *= binv;
	    element[1] *= binv;
	    element[2] *= binv;
	    return *this;
    }

    // Input / Output
    friend std::ostream & operator << (std::ostream & , const vector & );
	  friend std::istream & operator >> (std::istream & , vector & );
};

inline std::ostream & operator << ( std::ostream & s, const vector & v )
{
  return s << v.element[0] << "  "
           << v.element[1] << "  "
           << v.element[2];
}

inline std::istream & operator >> ( std::istream & s, vector & v )
{
  s >> v.element[0] >> v.element[1] >> v.element[2];
  return s;
}

inline real square ( vector v ) {return v * v;}
inline real abs ( vector v )    {return sqrt(v * v);}

inline vector operator + ( real b, const vector & v )
{
  return vector(b + v.element[0],
			          b + v.element[1],
			          b + v.element[2]);
}

inline vector operator + ( const vector & v, real b )
{
  return vector(b + v.element[0],
	              b + v.element[1],
			          b + v.element[2]);
}
inline vector operator * ( real b, const vector & v )
{
  return vector(b * v.element[0],
		            b * v.element[1],
			          b * v.element[2]);
}
inline vector operator * ( const vector & v, real b )
{
  return vector(b * v.element[0],
			          b * v.element[1],
			          b * v.element[2]);
}
inline vector operator / ( const vector & v, real b )
{
  return vector(v.element[0] / b,
	              v.element[1] / b,
			          v.element[2] / b);
}
#endif
