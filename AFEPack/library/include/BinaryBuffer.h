/**
 * @file   BinaryBuffer.h
 * @author Ruo Li <rli@aztec>
 * @date   Thu Sep 24 09:46:48 2009
 * 
 * @brief  二进制缓冲区及其上构造的流。
 * 
 * 
 */

#ifndef __BinaryBuffer_h__
#define __BinaryBuffer_h__

#include <vector>
#include "Geometry.h"

template <typename CHAR = char,
  template <typename T, typename ALLOC = std::allocator<T> > class VECTOR = std::vector>
  class BinaryBuffer : public VECTOR<CHAR> {
 public:
 typedef CHAR char_t;
 typedef VECTOR<CHAR> _Base;
 typedef typename _Base::iterator iterator;
 typedef typename _Base::const_iterator const_iterator;

 BinaryBuffer() : _Base() {}
 ~BinaryBuffer() {}

 /** 
  * 返回当前缓存的起始地址，供MPI传递函数使用
  */
 void * start_address() const { return (void *)&(*this)[0]; }

};

namespace AFEPack {

  template <typename BUFFER>
    class stream_base {
  protected:
    typedef typename BUFFER::char_t char_t;
    typedef typename BUFFER::const_iterator const_iterator;
  public:
    /// 空构造函数
  stream_base() : _buf(NULL) {}
    /// 提供一个数据缓存的构造函数 
  stream_base(BUFFER& buf) : _buf(&buf) {}
    /// 析构函数
    ~stream_base() {}

  public:
    /** 
     * 提供数据缓存对象
     * 
     * @param buf 
     */
    void set_buffer(BUFFER& buf) { _buf = &buf; }

    /** 
     * 查询得到数据缓存对象
     * @return  数据缓存
     */
    const BUFFER& get_buffer() const { return *_buf; }

  protected:
    /// 被操作的缓存对象
    BUFFER * _buf;
  };

  /// 支持数据写入缓存。
  template <typename BUFFER = BinaryBuffer<> >
    class ostream : public stream_base<BUFFER>{
  private:
  typedef stream_base<BUFFER> _Base;
  typedef typename _Base::char_t char_t;
  using _Base::_buf;
  public:
  ostream() : _Base() {}
  ostream(BUFFER& buf) : _Base(buf) {}
  ~ostream() {}

  /** 
   * 将类型为T数据t 放到当前缓存向量的末尾注意这里的T只支持基本数据类型的写
   * 入，其他任何非c++基本类型的写入需要用户自己给定写入格式。
   * 
   * @param t 待写入缓存的数据
   */
  template <class T>
  void encode(const T& t) {
    int n = sizeof(t);
    const char_t * ptr = (const char_t *) (&t);
    _buf->insert(_buf->end(), ptr, ptr + n);
  }

  void encode_binary(void * data, int n) {
    const char_t * ptr = (const char_t *)data;
    _buf->insert(_buf->end(), ptr, ptr + n); 
  }
  };

  /// 支持数据从缓存中读取。此对象构造完成后，目前只能从头到尾遍历读取一次。若要多次读取，需要重新构造流对象。
  template <typename BUFFER = BinaryBuffer<> >
    class istream : public stream_base<const BUFFER>{
  private:
  typedef stream_base<const BUFFER> _Base;
  typedef typename _Base::char_t char_t;
  typedef typename _Base::const_iterator iterator;
  using _Base::_buf;
  public:
  istream() : _Base() {}
  istream(const BUFFER& buf) : _Base(buf) { reset(); }
  ~istream() {}

  /// 设定作为数据源的缓存对象
  void set_buffer(const BUFFER& buf) {
    _Base::set_buffer(buf);
    reset();
  }

  /** 
   *  从当前向量读取数据到变量t中。这里同样只支持c++基本数据类型，若有需要处
   *  理其他类型的数据，用户需要提供带有特定数据类型参数的同名函数。
   * 
   * @param t 被读取的数据，同时被转换成类型T
   */
  template <class T>
  void decode(T& t) {
    std::size_t n = sizeof(T);
    std::copy(_pos, _pos + n, (char_t *)(&t));
    _pos += n;
  }
  void decode_binary(void * data, int n) {
    std::copy(_pos, _pos + n, (char_t *)data);
    _pos += n;
  }

  private:
  /** 
   * 将指示子放到数据缓存的开始位置，准备读取缓存数据
   * 
   */
  void reset() {
    if (_buf != NULL) {
      _pos = _buf->begin();
    }
  }

  private:
  /// 读取缓存是逐步进行的，因此需要一个迭代子指示当前的位置
  iterator _pos;
  }; 
  
  template <class OSTREAM, class T>
    void encode(OSTREAM& os, const T& t) {
    os.encode(t);
  }
  template <class ISTREAM, class T>
    void decode(ISTREAM& is, T& t) {
    is.decode(t);
  }

  template <class OSTREAM, class T>
    void encode_binary(OSTREAM& os, const T& t) {
    os.encode_binary(t);
  }
  template <class ISTREAM, class T>
    void decode_binary(ISTREAM& is, T& t) {
    is.decode_binary(t);
  }

  template <class OSTREAM, class T>
    OSTREAM& operator<<(OSTREAM& os, const T& t){
    os.encode(t);
    return os;
  }
  template <class ISTREAM, class T>
    ISTREAM& operator>>(ISTREAM& is, T& t){
    is.decode(t);
    return is;
  }



  template <class BUF, class T>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<T>& t) {
    u_int n = t.size();
    os << n;
    for (u_int i = 0;i < n;++ i) {
      os << t[i];
    }
    return os;
  }
  template <class BUF, class T>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<T>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    for (u_int i = 0;i < n;++ i) {
      is >> t[i];
    }
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<double>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t[0]), sizeof(double)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<double>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    is.decode_binary(&t[0], sizeof(double)*n);
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<int>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t[0]), sizeof(int)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<int>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    is.decode_binary(&t[0], sizeof(int)*n);
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<u_int>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t[0]), sizeof(u_int)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<u_int>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    is.decode_binary(&t[0], sizeof(u_int)*n);
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const Vector<double>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t(0)), sizeof(double)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, Vector<double>& t) {
    u_int n;
    is >> n;
    t.reinit(n);
    is.decode_binary(&t(0), sizeof(double)*n);
    return is;
  }

  template <class BUF, int DIM>
    ostream<BUF>& operator<<(ostream<BUF>& os, const Point<DIM>& t) {
    os.encode_binary((void *)(&t[0]), sizeof(double)*DIM);
    return os;
  }
  template <class BUF, int DIM>
    istream<BUF>& operator>>(istream<BUF>& is, Point<DIM>& t) {
    is.decode_binary(&t[0], sizeof(double)*DIM);
    return is;
  }

  template <class BUF, typename CHAR, template <typename T = char, 
    typename ALLOC = std::allocator<T> > class VEC>
    ostream<BUF>& operator<<(ostream<BUF>& os, const BinaryBuffer<CHAR,VEC>& buf) {
    const VEC<CHAR>& v(buf);
    os << v;
    return os;
  }
  template <class BUF, typename CHAR, template <typename T = char, 
    typename ALLOC = std::allocator<T> > class VEC>
    istream<BUF>& operator>>(istream<BUF>& is, BinaryBuffer<CHAR,VEC>& buf) {
    VEC<CHAR>& v(buf);
    is >> v;
    return is;
  }
}

#endif // __BinaryBuffer_h__

/**
 * End of File
 * 
 */
