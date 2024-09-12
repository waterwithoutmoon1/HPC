#include "atoms.h"
#include "berendsen_thermostat.h"
#include <gtest/gtest.h>

TEST(BerendsenThermostatTest, NoTemperatureChange) {
    // Create atoms with known velocities
    const size_t n = 5;

    Positions_t positions(3, n);
    positions << 0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0;

    Velocities_t velocities(3, n);
    velocities.setOnes();

    Atoms atoms{positions, velocities};

    // Set up the test environment
    double target_temperature = 1.0;
    double timestep = 0.001;
    double relaxation_time = 0.1;
    double kB = 1;
    double mass = 1;

    double temperature = (2.0 * compute_kinetic_energy(atoms.velocities, mass)) / (3.0 *  atoms.nb_atoms() * kB);
    // Check the temperature before applying the thermostat
    EXPECT_NEAR(temperature, target_temperature, 1e-6);

    // Apply the Berendsen thermostat
    apply_berendsen_thermostat(atoms, kB, mass, target_temperature, relaxation_time,timestep);

    // Check that velocities are scaled correctly
    EXPECT_NEAR(atoms.velocities(0, 0), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 0), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 0), 1, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 1), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 1), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 1), 1, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 2), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 2), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 2), 1, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 3), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 3), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 3), 1, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 4), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 4), 1, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 4), 1, 1e-6);
}

TEST(BerendsenThermostatTest, TemperatureDecrease) {
    // Create atoms with known velocities
    const size_t n = 5;

    Positions_t positions(3, n);
    positions << 0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0;

    Velocities_t velocities(3, n);
    velocities.setOnes();

    Atoms atoms{positions, velocities};

    // Set up the test environment
    double target_temperature = 0.9;
    double timestep = 0.001;
    double relaxation_time = 0.1;
    double kB = 1;
    double mass = 1;

    double temperature = (2.0 * compute_kinetic_energy(atoms.velocities, mass)) / (3.0 *  atoms.nb_atoms() * kB);
    // Check the temperature before applying the thermostat
    EXPECT_NEAR(temperature, target_temperature, 0.1);

    // Apply the Berendsen thermostat
    apply_berendsen_thermostat(atoms, kB, mass, target_temperature, relaxation_time,timestep);

    // Check that velocities are scaled correctly
    EXPECT_NEAR(atoms.velocities(0, 0), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 0), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 0), 0.9994998749374, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 1), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 1), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 1), 0.9994998749374, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 2), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 2), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 2), 0.9994998749374, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 3), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 3), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 3), 0.9994998749374, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 4), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 4), 0.9994998749374, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 4), 0.9994998749374, 1e-6);
}

TEST(BerendsenThermostatTest, TemperatureIncrease) {
    // Create atoms with known velocities
    const size_t n = 5;

    Positions_t positions(3, n);
    positions << 0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0;

    Velocities_t velocities(3, n);
    velocities.setOnes();

    Atoms atoms{positions, velocities};

    // Set up the test environment
    double target_temperature = 10.0;
    double timestep = 0.001;
    double relaxation_time = 0.1;
    double kB = 1;
    double mass = 1;

    double temperature = (2.0 * compute_kinetic_energy(atoms.velocities, mass)) / (3.0 *  atoms.nb_atoms() * kB);
    // Check the temperature before applying the thermostat
    EXPECT_NEAR(temperature, target_temperature, 9);

    // Apply the Berendsen thermostat
    apply_berendsen_thermostat(atoms, kB, mass, target_temperature, relaxation_time,timestep);

    // Check that velocities are scaled correctly
    EXPECT_NEAR(atoms.velocities(0, 0), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 0), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 0), 1.044030650891, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 1), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 1), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 1), 1.044030650891, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 2), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 2), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 2), 1.044030650891, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 3), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 3), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 3), 1.044030650891, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 4), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 4), 1.044030650891, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 4), 1.044030650891, 1e-6);
}