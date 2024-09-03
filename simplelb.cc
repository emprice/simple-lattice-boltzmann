#include "array2d.hh"
#include "array3d.hh"

#include <omp.h>
#include <hdf5.h>
#include <random>

struct SimpleLatticeBoltzmann
{
    SimpleLatticeBoltzmann(int nx, int ny);

    void set_initial_condition();
    void fill_moments();
    void step();

    Array2d<double> const& rho() const;
    Array2d<double> const& rho_ux() const;
    Array2d<double> const& rho_uy() const;

    protected:
        int m_nx;    ///< Size of mesh in x
        int m_ny;    ///< Size of mesh in y

        Array2d<double> m_rho;      ///< Current density
        Array2d<double> m_rhoux;    ///< Current x-momentum
        Array2d<double> m_rhouy;    ///< Current y-momentum

        Array3d<double> m_fcur;     ///< Current distribution
        Array3d<double> m_feq;      ///< Equilibrium distribution
        Array3d<double> m_fstar;    ///< Post-collision equilibrium distribution

        static constexpr int npop = 9;   ///< Number of particle populations

        /// Population velocity in x-coordinate
        static constexpr int cx[] = { 0, +1, 0, -1, 0, +1, -1, -1, +1 };

        /// Population velocity in y-coordinate
        static constexpr int cy[] = { 0, 0, +1, 0, -1, +1, +1, -1, -1 };

        /// Population weights, expressed as fraction of 36
        static constexpr int w[] = { 16, 4, 4, 4, 4, 1, 1, 1, 1 };
};

SimpleLatticeBoltzmann::SimpleLatticeBoltzmann(int nx, int ny) :
    m_nx(nx), m_ny(ny),
    m_rho(m_nx, m_ny), m_rhoux(m_nx, m_ny), m_rhouy(m_nx, m_ny),
    m_fcur(m_nx, m_ny, npop), m_feq(m_nx, m_ny, npop), m_fstar(m_nx, m_ny, npop)
{
    // Set the initial population values
    set_initial_condition();

    // Ensure the density and momenta are valid even if no steps are taken
    fill_moments();
}

/**
 * @brief Fill a randomly-perturbed initial condition (for demo purposes)
 */
void SimpleLatticeBoltzmann::set_initial_condition()
{
    std::normal_distribution<double> normal(0, 0.1);
    std::uniform_real_distribution<double> uniform(0, 0.1);

    // Zero the initial distribution function; we can collapse these loops
    // for better performance
#pragma omp parallel for collapse(3)
    for (int i = 0; i < m_nx; ++i)
    {
        for (int j = 0; j < m_ny; ++j)
        {
            for (int k = 0; k < npop; ++k)
            {
                m_fcur(i,j,k) = 0;
            }
        }
    }

    // Set the non-zero population with a perturbation, being careful to
    // avoid repeat values from non-threadsafe generators
#pragma omp parallel
    {
        // A simple approach is to seed the generator with the thread number,
        // but we can only do it inside a parallel region, before the for loop
        std::default_random_engine generator(omp_get_thread_num());

#pragma omp for collapse(2)
        for (int i = 0; i < m_nx; ++i)
        {
            for (int j = 0; j < m_ny; ++j)
            {
                m_fcur(i,j,0) = 1 + normal(generator);
                m_fcur(i,j,1) = 2 + uniform(generator);
            }
        }
    }
}

/**
 * @brief Fill the density, x-momentum, and y-momentum for the current
 * distribution stored in `m_fcur`
 */
void SimpleLatticeBoltzmann::fill_moments()
{
    // Compute density and momenta from current populations
#pragma omp parallel for collapse(2)
    for (int i = 0; i < m_nx; ++i)
    {
        for (int j = 0; j < m_ny; ++j)
        {
            double rho_val = 0, rhoux_val = 0, rhouy_val = 0;

            // It would be nice to collapse this loop, but it's not perfectly
            // nested; with more populations, we could potentially parallelize
            // the inner loop, but this is fine for now.
            for (int k = 0; k < npop; ++k)
            {
                // Only read from memory once (more efficient than reading the
                // same value three times)
                double f_val = m_fcur(i,j,k);

                rho_val  += f_val;
                rhoux_val += f_val * cx[k];
                rhouy_val += f_val * cy[k];
            }

            // Only write to memory once (more efficient writing in the loop)
            m_rho(i,j) = rho_val;
            m_rhoux(i,j) = rhoux_val;
            m_rhouy(i,j) = rhouy_val;
        }
    }
}

void SimpleLatticeBoltzmann::step()
{
    // In a more complex implementation, these should be specified by the
    // user; for now, these constant values suffice
    static constexpr double dt = 1;
    static constexpr double tau = 3;
    static const double dt_over_tau = dt / tau;

    // Fill the equilibrium populations; some of the computations are repeated,
    // but this allows for tighter nesting of the for loops so we can collapse
    // all levels
#pragma omp parallel for collapse(3)
    for (int i = 0; i < m_nx; ++i)
    {
        for (int j = 0; j < m_ny; ++j)
        {
            for (int k = 0; k < npop; ++k)
            {
                // Load from memory once
                double rho_val = m_rho(i,j);

                // Precompute these values for efficiency
                double ux = m_rhoux(i,j) / rho_val;
                double uy = m_rhouy(i,j) / rho_val;

                // Compute u.u (velocity dotted with itself)
                double udotu = ux * ux + uy * uy;

                // Convert the integer weight to its actual value
                double weight = w[k] / 36.;

                // Compute u.c (velocity dotted with soundspeed)
                double udotc = ux * cx[k] + uy * cy[k];

                // Assuming the soundspeed constant is unity
                m_feq(i,j,k) = rho_val * weight *
                    (1 + udotc + 0.5 * udotc * udotc - 0.5 * udotu);
            }
        }
    }

    // Trivially parallelized vector operation converts current distribution
    // and equilibrium distribution to post-collision distribution
#pragma omp parallel for collapse(3)
    for (int i = 0; i < m_nx; ++i)
    {
        for (int j = 0; j < m_ny; ++j)
        {
            for (int k = 0; k < npop; ++k)
            {
                m_fstar(i,j,k) = m_fcur(i,j,k) * (1 - dt_over_tau) +
                    m_feq(i,j,k) * dt_over_tau;
            }
        }
    }

    // Thanks to above precomputation, the streaming step is a simple periodic
    // shift in indices
#pragma omp parallel for collapse(3)
    for (int i = 0; i < m_nx; ++i)
    {
        for (int j = 0; j < m_ny; ++j)
        {
            for (int k = 0; k < npop; ++k)
            {
                // Source indices need to be offset and made periodic
                int ii = (i - cx[k] + m_nx) % m_nx;
                int jj = (j - cy[k] + m_ny) % m_ny;

                // New "current" distribution is a shifted post-collision
                // distribution
                m_fcur(i,j,k) = m_fstar(ii,jj,k);
            }
        }
    }

    // Update the moments now that the current distribution has been updated
    fill_moments();
}

/// Getter for current density array
Array2d<double> const& SimpleLatticeBoltzmann::rho() const { return m_rho; }

/// Getter for current x-momentum array
Array2d<double> const& SimpleLatticeBoltzmann::rho_ux() const { return m_rhoux; }

/// Getter for current y-momentum array
Array2d<double> const& SimpleLatticeBoltzmann::rho_uy() const { return m_rhouy; }

/**
 * @brief Write the given array to an open dataset
 * @param[in] arr Array to write
 * @param[in] dset Open dataset
 * @param[in] mspace Memory space descriptor
 * @param[in] dims Desired dimensions of the dataset after writing
 * @param[in] offset Offset in the file to write to
 * @param[in] chunk Size of the data to write
 */
void write_data(Array2d<double> const& arr, hid_t dset,
    hid_t mspace, hsize_t dims[], hsize_t offset[], hsize_t chunk[])
{
    H5Dextend(dset, dims);
    hid_t fspace = H5Dget_space(dset);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, chunk, nullptr);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT,
        const_cast<double*>(arr.data()));
    H5Sclose(fspace);
}

int main(int argc, char *argv[])
{
    static constexpr int nx = 100, ny = 100;

    SimpleLatticeBoltzmann lbm(nx, ny);

    // Open the output file
    hid_t fid = H5Fcreate("simple_lbm.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t cur_dims[3] = { 0, nx, ny };
    hsize_t chunk_dims[3] = { 1, nx, ny };
    hsize_t max_dims[3] = { H5S_UNLIMITED, nx, ny };
    hsize_t offset[3] = { 0, 0, 0 };

    // Set up for an unlimited-dimension write; useful for codes where number
    // of steps is controlled by the user
    hid_t spid = H5Screate_simple(3, cur_dims, max_dims);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 3, chunk_dims);

    // Using native type in the file means no conversion needs to be done
    // on-the-fly, which is faster
    hid_t dset_rho = H5Dcreate(fid, "rho", H5T_NATIVE_DOUBLE,
        spid, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    hid_t dset_rho_ux = H5Dcreate(fid, "rho_ux", H5T_NATIVE_DOUBLE,
        spid, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    hid_t dset_rho_uy = H5Dcreate(fid, "rho_uy", H5T_NATIVE_DOUBLE,
        spid, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    H5Pclose(dcpl);

    // In memory, we only store one 2d slab
    hid_t mspace = H5Screate_simple(3, chunk_dims, nullptr);
    H5Sselect_all(mspace);

    // Write the initial condition
    cur_dims[0] += 1;
    write_data(lbm.rho(), dset_rho, mspace, cur_dims, offset, chunk_dims);
    write_data(lbm.rho_ux(), dset_rho_ux, mspace, cur_dims, offset, chunk_dims);
    write_data(lbm.rho_uy(), dset_rho_uy, mspace, cur_dims, offset, chunk_dims);
    offset[0] += 1;

    for (int step = 0; step < 100; ++step)
    {
        // Take a single simulation step
        lbm.step();

        // Write the result of the last step
        cur_dims[0] += 1;
        write_data(lbm.rho(), dset_rho, mspace, cur_dims, offset, chunk_dims);
        write_data(lbm.rho_ux(), dset_rho_ux, mspace, cur_dims, offset, chunk_dims);
        write_data(lbm.rho_uy(), dset_rho_uy, mspace, cur_dims, offset, chunk_dims);
        offset[0] += 1;
    }

    // Clean up file resources
    H5Sclose(mspace);
    H5Dclose(dset_rho);
    H5Dclose(dset_rho_ux);
    H5Dclose(dset_rho_uy);
    H5Fclose(fid);

    return 0;
}

// vim: set ft=cpp.doxygen:
