/**
 * @file   Migration.h
 * @author Ruo Li <rli@aztec>
 * @date   Thu Sep 24 11:51:29 2009
 * 
 * @brief  数据迁移
 * 
 * 
 */

#ifndef __Migration_h__
#define __Migration_h__

#include <map>
#include "BinaryBuffer.h"
#include "TemplateElement.h"

/**
 * 在进行整个环境的序列化以前，我们先构造序列化的环境。我们设定基本前提
 * 为：假定进程上仅仅有一棵几何遗传树，在此树上多个非正则网格，在每个非
 * 正则网格上建立了一个正则网格，在正则网格上则建立了多个有限元空间，每
 * 个有限元空间上有多个有限元函数。
 *
 * 由于在进行网格的拆分和合并的时候，正则网格需要重建，从而有限元空间也
 * 需要进行重建，为了能够对有限元函数的信息进行迁移，我们在猜分以前，将
 * 有限元函数的信息存储到相应的几何遗传树的几何体上去。几何遗传树中的几
 * 何体继承了 HGeomeryBuffer 这个类，从而我们可以将信息存储在这个类的缓
 * 冲区中。比如从有限元函数的出发，每个自由度都能够找到自身所依附的正则
 * 网格中的几何体，而每个正则网格中的几何体则需要能够找寻到在正则化的过
 * 程中相应的几何遗传树中的几何体，我们对正则网格进行了改进，使得这个检
 * 索成为可能。这个是通过在正则网格中加入了h_geometry_ptr数组来完成的。
 *
 * 在此基础上，我们假设在每一个几何体上附着了一组数据，这些数据按照数
 * 据 ID 被分成组，比如一个有限元函数的数据我们可以分为一组。一个几何
 * 体上，属于同一个数据 ID 的数据则必须由用户自己负责分析，也就是说，
 * 这组数据本身必须在上下文中是可以自解释的。具体参考 export_fe_func
 * 和 import_fe_func 这对函数的实现。
 *
 */
namespace Migration {

  typedef int data_id_t; /// 数据 ID 的类型：为有符号整型
  typedef std::string data_name_t; /// 数据名称类型：字符串
  typedef BinaryBuffer<> data_buffer_t; /// 数据缓冲区类型
  /**
   * 数据 ID 到缓冲区映射表，这个是用在最终存储数据位置的类型。
   */
  typedef std::map<data_id_t, data_buffer_t> buffer_t; 

  data_id_t name_to_id(const data_name_t& dn); /// 从数据名称到ID的转换
  data_id_t register_data_name(const data_name_t& dn, bool); /// 登记数据名称
  void initialize(); /// 初始化数据迁移环境
  bool is_valid(const data_id_t&); /// 判断一个数据ID是否合法

  namespace details {
    /**
     * 获取树结构中几何体上的数据缓冲区。对于存储的情形，如果相应的缓冲
     * 区不存在，则新建立改缓冲区。
     * 
     * @param geo 树结构中的几何体 
     * @param is_save 存储还是读入
     * 
     * @return 相应的缓冲区
     */
    template <class HGEO> BinaryBuffer<>& 
      get_buffer(HGEO& geo,
                 const data_id_t& data_id,
                 bool is_save) {
      buffer_t& buffer = geo.buffer;
      if (is_save) {
        return buffer[data_id];
      } else {
        buffer_t::iterator it = buffer.find(data_id);
        if (it == buffer.end()) {
          return *((BinaryBuffer<> *)NULL);
        }
        return it->second;
      }
    }

    /**
     * 获取几何体上的数据缓冲区。对于存储的情形，如果相应的缓冲区不存在，
     * 则新建立改缓冲区。
     * 
     * @param mesh 正则网格
     * @param data_id 数据 ID
     * @param dimension 几何体维数
     * @param geo_idx 几何体序号
     * @param is_save 存储还是读入
     * 
     * @return 相应的缓冲区
     */
    template <class MESH> BinaryBuffer<>& 
      get_buffer(MESH& mesh,
                 const data_id_t& data_id,
                 u_int dimension,
                 u_int geo_idx,
                 bool is_save) {
      buffer_t& buffer = mesh.h_geometry(dimension, geo_idx)->buffer;
      if (is_save) {
        return buffer[data_id];
      } else {
        buffer_t::iterator it = buffer.find(data_id);
        if (it == buffer.end()) {
          return *((BinaryBuffer<> *)NULL);
        }
        return it->second;
      }
    }

  }

  /** 
   * 获取树结构中的几何体上的一个输出流，用来输出数据。
   * 
   * @param geo 树结构中的几何体
   * @param data_id 数据 ID
   * @param os 获得的输出流
   */
  template <class HGEO, class STREAM> void
    get_export_stream(HGEO& geo,
                      const data_id_t& data_id,
                      STREAM& os) {
    os.set_buffer(details::get_buffer(geo, data_id, true));
  }
  
  /** 
   * 获取树结构中几何体上的一个输入流，用来输出数据。
   * 
   * @param geo 树结构中的几何体
   * @param data_id 数据 ID
   * @param is 获得的输入流
   */
  template <class HGEO, class STREAM> void
    get_import_stream(HGEO& geo,
                      const data_id_t& data_id,
                      STREAM& is) {
    is.set_buffer(details::get_buffer(geo, data_id, false));
  }

  /** 
   * 获取网格中的几何体上的一个输出流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param data_id 数据 ID
   * @param dimension 几何体维数
   * @param geo_idx 几何体序号
   * @param os 获得的输出流
   */
  template <class MESH, class STREAM> void
    get_export_stream(MESH& mesh,
                      const data_id_t& data_id,
                      u_int dimension,
                      u_int geo_idx,
                      STREAM& os) {
    os.set_buffer(details::get_buffer(mesh, data_id, dimension, geo_idx, true));
  }
  
  /** 
   * 获取网格中几何体上的一个输入流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param data_id 数据 ID
   * @param dimension 几何体维数
   * @param geo_idx 几何体序号
   * @param is 获得的输入流
   */
  template <class MESH, class STREAM> void
    get_import_stream(MESH& mesh,
                      const data_id_t& data_id,
                      u_int dimension,
                      u_int geo_idx,
                      STREAM& is) {
    is.set_buffer(details::get_buffer(mesh, data_id, dimension, geo_idx, false));
  }

  /** 
   * 获取有限元空间上的一个自由度的输出流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param sp 有限元空间
   * @param data_id 数据 ID
   * @param dof 自由度编号
   * @param os 获得的输出流
   */
  template <class MESH, class SP, class STREAM> void
    get_dof_export_stream(MESH& mesh,
                          SP& sp,
                          const data_id_t& data_id,
                          int dof,
                          STREAM& os) {
    const DOFIndex& di = sp.dofIndex(dof);
    get_export_stream(mesh, data_id, di.dimension, di.geometry_index, os);
  }

  /** 
   * 获取有限元空间上的一个自由度的输入流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param sp 有限元空间
   * @param data_id 数据 ID
   * @param dof 自由度编号
   * @param is 获得的输入流
   */
  template <class MESH, class SP, class STREAM> void
    get_dof_import_stream(MESH& mesh,
                          SP& sp,
                          const data_id_t& data_id,
                          int dof,
                          STREAM& is) {
    const DOFIndex& di = sp.dofIndex(dof);
    get_import_stream(mesh, data_id, di.dimension, di.geometry_index, is);
  }

} // end of namespace Migration

#endif // __Migration_h__

/**
 * end of file
 * 
 */
 
