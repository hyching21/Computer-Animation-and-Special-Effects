#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first because it is very simple. Then you can check whether your collision is correct or
    //   not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Review "ODE_basics.pptx" from p.15 - p.16

    float dt = particleSystem.deltaTime;

    // update elevator
    Eigen::Vector3f elevatorPos = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorVel = particleSystem.elevatorTerrain->getVelocity();
    Eigen::Vector3f elevatorAcc = particleSystem.elevatorTerrain->getAcceleration();

    Eigen::Vector3f newElevatorVel = elevatorVel + elevatorAcc * dt;
    Eigen::Vector3f newElevatorPos = elevatorPos + elevatorVel * dt;

    particleSystem.elevatorTerrain->setVelocity(newElevatorVel);
    particleSystem.elevatorTerrain->setPosition(newElevatorPos);
    particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());

    // update jelly
    for(int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);

        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            // x = x + v * dt
            particle.addPosition(particle.getVelocity() * particleSystem.deltaTime);
            // v = v + a * dt
            particle.addVelocity(particle.getAcceleration() * particleSystem.deltaTime);
            // Clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_implicit.pptx" from p.18 - p.19

    Eigen::Vector3f elevatorPosBackup = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorVelBackup = particleSystem.elevatorTerrain->getVelocity();
    int elevatorCounterBackup = particleSystem.elevatorCounter;

    for(int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        Jelly backupJelly = *jelly;        
        float dt = particleSystem.deltaTime;

        // reset force
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            particle.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());

        particleSystem.computeElevatorForce();
        particleSystem.computeJellyForce(*jelly);
        particleSystem.elevatorCounter = elevatorCounterBackup;

        // update elevator
        Eigen::Vector3f elevatorVel = particleSystem.elevatorTerrain->getVelocity();
        Eigen::Vector3f elevatorAcc = particleSystem.elevatorTerrain->getAcceleration();
        Eigen::Vector3f newElevatorVel = elevatorVelBackup + elevatorAcc * dt;
        Eigen::Vector3f newElevatorPos = elevatorPosBackup + newElevatorVel * dt;

        particleSystem.elevatorTerrain->setPosition(newElevatorPos);
        particleSystem.elevatorTerrain->setVelocity(newElevatorVel);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());

        // update jelly
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Particle& backupParticle = backupJelly.getParticle(j);

            Eigen::Vector3f vel = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();

            // v(n+1) = v(n) + F(n+1)/m * dt
            Eigen::Vector3f newVel = backupParticle.getVelocity() + acc * dt;
            // x(n+1) = x(n) + v(n+1) * dt
            Eigen::Vector3f newPos = backupParticle.getPosition() + newVel * dt;
            particle.setPosition(newPos);
            particle.setVelocity(newVel);
            // Clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_basics.pptx" from p.18 - p.19

    Eigen::Vector3f elevatorPosBackup = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorVelBackup = particleSystem.elevatorTerrain->getVelocity();
    int elevatorCounterBackup = particleSystem.elevatorCounter;

    for(int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        Jelly backupJelly = *jelly;        
        float dt = particleSystem.deltaTime;

        // reset force
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f vel = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();
            Eigen::Vector3f midpos = particle.getPosition() + 0.5 * vel * dt;
            Eigen::Vector3f midvel = particle.getVelocity() + 0.5 * acc * dt;
            particle.setPosition(midpos);
            particle.setVelocity(midvel);
            particle.setForce(Eigen::Vector3f::Zero());
        }
        
        Eigen::Vector3f ele_midpos = particleSystem.elevatorTerrain->getPosition() + 0.5 * particleSystem.elevatorTerrain->getVelocity() * dt;
        Eigen::Vector3f ele_midvel = particleSystem.elevatorTerrain->getVelocity() + 0.5 * particleSystem.elevatorTerrain->getAcceleration() * dt;
        particleSystem.elevatorTerrain->setPosition(ele_midpos);
        particleSystem.elevatorTerrain->setVelocity(ele_midvel);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());

        particleSystem.computeElevatorForce();
        particleSystem.computeJellyForce(*jelly);
        particleSystem.elevatorCounter = elevatorCounterBackup;

        // update elevator
        Eigen::Vector3f elevatorVel = particleSystem.elevatorTerrain->getVelocity();
        Eigen::Vector3f elevatorAcc = particleSystem.elevatorTerrain->getAcceleration();
        Eigen::Vector3f newElevatorVel = elevatorVelBackup + elevatorAcc * dt;
        Eigen::Vector3f newElevatorPos = elevatorPosBackup + elevatorVel * dt;

        particleSystem.elevatorTerrain->setPosition(newElevatorPos);
        particleSystem.elevatorTerrain->setVelocity(newElevatorVel);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());

        // update jelly
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Particle& backupParticle = backupJelly.getParticle(j);

            Eigen::Vector3f vel = particle.getVelocity();
            Eigen::Vector3f acc = particle.getAcceleration();
            // v(n+1) = v(n) + a(n+0.5) * dt
            Eigen::Vector3f newVel = backupParticle.getVelocity() + acc * dt;
            // x(n+1) = x(n) + v(n+0.5) * dt
            Eigen::Vector3f newPos = backupParticle.getPosition() + vel * dt;
            particle.setPosition(newPos);
            particle.setVelocity(newVel);
            // Clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    float dt = particleSystem.deltaTime;
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. StateStep struct is just a hint, you can use whatever you want.
    //   6. Review "ODE_basics.pptx" from p.21

    // #### For Elevator ####
    Eigen::Vector3f elevatorPosBackup = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorVelBackup = particleSystem.elevatorTerrain->getVelocity();
    int elevatorCounterBackup = particleSystem.elevatorCounter;
    StateStep e_k1, e_k2, e_k3, e_k4, e_total;
    // k1
    particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    particleSystem.computeElevatorForce();
    e_k1.deltaVel = particleSystem.elevatorTerrain->getAcceleration() * dt;
    e_k1.deltaPos = elevatorVelBackup * dt;
    // prepare for k2 state
    particleSystem.elevatorTerrain->setVelocity(elevatorVelBackup + 0.5 * e_k1.deltaVel);
    particleSystem.elevatorTerrain->setPosition(elevatorPosBackup + 0.5 * e_k1.deltaPos);
    particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    particleSystem.elevatorCounter = elevatorCounterBackup;
    particleSystem.computeElevatorForce();

    // k2
    e_k2.deltaVel = particleSystem.elevatorTerrain->getAcceleration() * dt;
    e_k2.deltaPos = particleSystem.elevatorTerrain->getVelocity() * dt;
    // prepare for k3 state
    particleSystem.elevatorTerrain->setVelocity(elevatorVelBackup + 0.5 * e_k2.deltaVel);
    particleSystem.elevatorTerrain->setPosition(elevatorPosBackup + 0.5 * e_k2.deltaPos);
    particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    particleSystem.elevatorCounter = elevatorCounterBackup;
    particleSystem.computeElevatorForce();

    // k3
    e_k3.deltaVel = particleSystem.elevatorTerrain->getAcceleration() * dt;
    e_k3.deltaPos = particleSystem.elevatorTerrain->getVelocity() * dt;
    // prepare for k4 state
    particleSystem.elevatorTerrain->setVelocity(elevatorVelBackup + e_k3.deltaVel);
    particleSystem.elevatorTerrain->setPosition(elevatorPosBackup + e_k3.deltaPos);
    particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    particleSystem.elevatorCounter = elevatorCounterBackup;
    particleSystem.computeElevatorForce();

    // k4
    e_k4.deltaVel = particleSystem.elevatorTerrain->getAcceleration() * dt;
    e_k4.deltaPos = particleSystem.elevatorTerrain->getVelocity() * dt;
    e_total = {(e_k1.deltaVel + 2.0 * e_k2.deltaVel + 2.0 * e_k3.deltaVel + e_k4.deltaVel) / 6.0f,
                        (e_k1.deltaPos + 2.0 * e_k2.deltaPos + 2.0 * e_k3.deltaPos + e_k4.deltaPos) / 6.0f};
    particleSystem.elevatorTerrain->setPosition(elevatorPosBackup + e_total.deltaPos);
    particleSystem.elevatorTerrain->setVelocity(elevatorVelBackup + e_total.deltaVel);
    particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    particleSystem.elevatorCounter = elevatorCounterBackup;

    // #### FOR Jelly ####
    for(int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        Jelly backupJelly = *jelly;        
        std::vector<Eigen::Vector3f> initPos, initVel;
        std::vector<StateStep> k1, k2, k3, k4;

        // k1
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            initPos.push_back(particle.getPosition());
            initVel.push_back(particle.getVelocity());

            StateStep step;
            step.deltaVel = particle.getAcceleration() * dt;
            step.deltaPos = particle.getVelocity() * dt;
            k1.push_back(step);

            // prepare for k2 state
            particle.addPosition(0.5 * step.deltaPos);
            particle.addVelocity(0.5 * step.deltaVel);
            particle.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*jelly);

        // k2
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            StateStep step;
            step.deltaVel = particle.getAcceleration() * dt;
            step.deltaPos = particle.getVelocity() * dt;
            k2.push_back(step);

            // prepare for k3 state
            particle.setPosition(initPos[j] + 0.5 * step.deltaPos);
            particle.setVelocity(initVel[j] + 0.5 * step.deltaVel);
            particle.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*jelly);

        // k3
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            StateStep step;
            step.deltaVel = particle.getAcceleration() * dt;
            step.deltaPos = particle.getVelocity() * dt;
            k3.push_back(step);

            // prepare for k4 state
            particle.setPosition(initPos[j] + step.deltaPos);
            particle.setVelocity(initVel[j] + step.deltaVel);
            particle.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*jelly);

        // k4
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            StateStep step, total;
            step.deltaVel = particle.getAcceleration() * dt;
            step.deltaPos = particle.getVelocity() * dt;
            k4.push_back(step);

            total.deltaPos = (k1[j].deltaPos + 2.0 * k2[j].deltaPos + 2.0 * k3[j].deltaPos + k4[j].deltaPos) / 6.0f;
            total.deltaVel = (k1[j].deltaVel + 2.0 * k2[j].deltaVel + 2.0 * k3[j].deltaVel + k4[j].deltaVel) / 6.0f;
            
            particle.setPosition(initPos[j] + total.deltaPos);
            particle.setVelocity(initVel[j] + total.deltaVel);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}
}  // namespace simulation
