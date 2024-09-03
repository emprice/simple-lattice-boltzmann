#ifndef ARRAY3D_HH
#define ARRAY3D_HH

#include <cassert>
#include <algorithm>

/**
 * @brief Container class for a basic 3d array
 * @tparam T Datatype to be stored
 */
template <typename T>
struct Array3d
{
    /// Default constructor for an empty, zero-size array
    Array3d() : m_data(nullptr), m_size1(0),
        m_size2(0), m_size3(0), m_owndata(false) { }

    /// Constructs a new array and allocates its memory
    Array3d(size_t size1, size_t size2, size_t size3) :
        m_size1(size1), m_size2(size2), m_size3(size3), m_owndata(true)
    {
        m_data = new T[m_size1*m_size2*m_size3];
    }

    /// Basic destructor to clean up heap-allocated memory
    virtual ~Array3d()
    {
        if (m_owndata) delete[] m_data;
    }

    /// Disabled copy constructor
    Array3d(Array3d const& other) = delete;

    /// Disabled copy assignment
    Array3d& operator=(Array3d const& other) = delete;

    /// Move constructor (no copies)
    Array3d(Array3d&& other) : Array3d()
    {
        std::swap(m_data, other.m_data);
        std::swap(m_size1, other.m_size1);
        std::swap(m_size2, other.m_size2);
        std::swap(m_size3, other.m_size3);
        std::swap(m_owndata, other.m_owndata);
    }

    /// Setter for element (i,j,k) of the array
    T& operator()(size_t i, size_t j, size_t k)
    {
        assert((i < m_size1) && (j < m_size2) && (k < m_size3));
        return m_data[(i*m_size2+j)*m_size3+k];
    }

    /// Getter for element (i,j,k) of the array
    T operator()(size_t i, size_t j, size_t k) const
    {
        assert((i < m_size1) && (j < m_size2) && (k < m_size3));
        return m_data[(i*m_size2+j)*m_size3+k];
    }

    /// Getter for the raw data buffer
    T *data() { return m_data; }

    /// Const getter for the raw data buffer
    const T *data() const { return m_data; }

    protected:
        T *m_data;          ///< Raw data pointer
        size_t m_size1;     ///< Array size on dimension 1
        size_t m_size2;     ///< Array size on dimension 2
        size_t m_size3;     ///< Array size on dimension 3
        bool m_owndata;     ///< Flag for data ownership
};

#endif      // ARRAY3D_HH

// vim: set ft=cpp.doxygen:
