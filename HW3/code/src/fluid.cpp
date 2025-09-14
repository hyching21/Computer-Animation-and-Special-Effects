#include "fluid.h"
#include <iostream>
void Fluid::integrate()
{
	Eigen::Vector2f gravityVec(0, -gravity * dt);

	for (int i = 0; i < num_particles; ++i)
	{
		particle_vel[i] += gravityVec;
		particle_pos[i] += dt * particle_vel[i];
	}
}

void Fluid::handleCollisions(bool lmbDown)
{
	Eigen::Vector2f minVal(cell_dim + particle_radius, cell_dim + particle_radius),
		maxVal(viewport_w - cell_dim - particle_radius, viewport_h - cell_dim - particle_radius);
	
	float min_obstacle_d2 = (obstacle_r + particle_radius) * (obstacle_r + particle_radius);
	Eigen::Vector2f obstacle_pos(obstacle_x, obstacle_y), obstacle_vel(obstacle_vx, obstacle_vy);

	for (int i = 0; i < num_particles; ++i)
	{
		Eigen::Vector2f& pos = particle_pos[i];
		Eigen::Vector2f& vel = particle_vel[i];

		// collision with obstacle
		if (lmbDown)
		{
			float d2_to_obstacle = (pos - obstacle_pos).squaredNorm();
			if (d2_to_obstacle < min_obstacle_d2)
				vel = obstacle_vel;
		}

		// collision with edges
		Eigen::Vector2f clamped = pos.cwiseMax(minVal).cwiseMin(maxVal);

		if (clamped.x() != pos.x()) vel.x() = 0.0f;
		if (clamped.y() != pos.y()) vel.y() = 0.0f;

		pos = clamped;
	}
}

void Fluid::particleRelaxation()
{
	/*TODO*/
	/*
	This function handles inter-particle collisions, here's what you'll have to implement:

	1.	Load all particles indices into relaxation cells

	2.	Implement the iterative solver:

		At each iteration, go through all relaxation cells.
		For each particle in the cell, compare its position to every other particle within a 3x3 region.
		If two particles are colliding, i.e. |d| < 2 * r, shift them away from each other so the two particles are 2 * r units apart.
		[DEBUGGING HINT: Two particles may completely overlap with each other!]

		You can directly apply the displacement vectors right after their calculation.


	Some of the variables you might need, but not limited to:

	Relaxation Cells:
		relaxation_cell_col				// Number of relaxation cell columns
		relaxation_cell_rows			// Number of relaxation cell rows
		relaxation_cell_particle_ids	// Index of particles in each cell
		relaxation_cell_dim				// Dimension of each relaxation cell (same width and height)


	Particle Information:
		particle_pos					// Position (x, y) of each particle
		particle_radius	
		num_particles
	*/

	// TODO: Assign particles to cells
	for (int i = 0; i < relaxation_cell_rows; ++i) {
		for (int j = 0; j < relaxation_cell_cols; ++j) {
			relaxation_cell_particle_ids[i][j].clear();
		}
	}

	for(int i = 0; i < num_particles; ++i) {
		Eigen::Vector2f& pos = particle_pos[i];
		int cell_x = static_cast<int>(pos.x() / relaxation_cell_dim);
		int cell_y = static_cast<int>(pos.y() / relaxation_cell_dim);

		if (cell_x >= 0 && cell_x < relaxation_cell_cols && cell_y >= 0 && cell_y < relaxation_cell_rows)
			relaxation_cell_particle_ids[cell_y][cell_x].push_back(i);
	}
	
	
	for (int iter = 0; iter < iterations; ++iter)
	{
		for (int i = 0; i < relaxation_cell_rows; ++i) for (int j = 0; j < relaxation_cell_cols; ++j)
		{
			// TODO: Perform particle relaxation
			for(int x = -1; x <= 1; ++x) for(int y = -1; y <= 1; ++y)
			{
				int cell_row = i + y; // row
				int cell_col = j + x; // col
				
				if (cell_col >= 0 && cell_col < relaxation_cell_cols && cell_row >= 0 && cell_row < relaxation_cell_rows){
					const auto& neighbor_particles = relaxation_cell_particle_ids[cell_row][cell_col];

					for (int p1 : relaxation_cell_particle_ids[i][j]) {
						for (int p2 : relaxation_cell_particle_ids[cell_row][cell_col]) {
							if (p1 >= p2) continue;

							Eigen::Vector2f d = particle_pos[p1] - particle_pos[p2];
							float dist = std::max(d.norm(), 1e-6f); 
							if (dist < 2 * particle_radius) {
								Eigen::Vector2f direction = d / dist;
								Eigen::Vector2f correction = 0.5f * direction * (2 * particle_radius - dist);
								particle_pos[p1] += correction;
								particle_pos[p2] -= correction;
							}
						}
					}
				}
			}
		}
	}
}

void Fluid::transferVelocities(bool to_cell)
{
	/*TODO*/
	/*
	This function has a parameter to_cell that determines whether the particle velocities should be transfered to the cells 
	or in the other direction.

	1.	Some setup you need to do when to_cell is true:

		Reset all cell_velocities.
		Assign each cell their proper cell type. A cell should be assigned the FLUID type if it contains any particles and AIR if not.
		SOLID cells are present on initialization and should not be modified.

	2.	In main loop:

		For each of x and y components, calculate the 4 cells that should contribute to the bilinear interpolation of each particle.
		From this step, we can obtain 4 weights, w_1 ~ w_4. These weights will be used in both scenarios (to_cell = true or false).

		When transfering particle velocities to the cells, add w_n * v to each of the 4 cell velocity components. 

		We should also accumulate w_n of each sampling point in a separate buffer to normalize the cell velocities after 
		iterating through all the particles.

		When transfering cell velocities back to the particles, before doing the usual bilinear interpolation, we need to test if 
		the velocities are valid.
		
		Let's take a look at an example, where we calculate the x velocity component
	
				|		|
		-------------------------
				|		|			vx_1: Vx(i, j)
			F	•->	A	•->	A		vx_2: Vx(i + 1, j)
		   Vx_4	|   	| Vx_3		vx_3: Vx(i + 1, j + 1)
		-------------------------	vx_4: Vx(i, j + 1)
				|	  ⦾	|			
			F	•->	F <-•	A		F: Fluid Cells
		   Vx_1	|		| Vx_2		A: Air Cells
		-------------------------
				|		|

		In this case, only vx_3 is invalid since the value it carries is derived from 2 AIR cells.
		Note, since the Vy is sampled at the center of the horizontal edge, we should test the upper and lower cells in that case.

		After this we can finally obtain the interpolated velocity. In the provided example, since Vx_3 is invalid, the interpolated velocity
		will be:

					w_1 * vx_1 + w_2 * vx_2 + w_4 * vx_4
		 Vx_pic =	------------------------------------
							  w_1 + w_2 + w_4

		We're almost done!
		The particle velocity was calculated using the PIC method. To obtain FLIP velocity use this formula:

							w_1 * (vx_1 - prev_v1) + w_2 * (vx_2 - prev_v2) + w_4 * (vx_4 - prev_v4)
		 Vx_flip =	Vx_p +	-----------------------------------------------------------------------
													   w_1 + w_2 + w_4

		We can finally set the particle velocity using a blend of PIC and FLIP velocities:

		 Vx_p = (1 - flip_ratio) * Vx_pic + flip_ratio * Vx_flip

		Don't forget to do this in the Y direction as well!

	3.	After the loop, if we're transfering cell velocities to the particle, we'll have to normalize each velocity component using 
		the stored normalizing values mentioned before.

		Backup cell velocities in prev_cell_velocities for FLIP calculation.


	Some of the variables you might need, but not limited to:

	flip_ratio					# used during pic/flip interpolation

	Particle Information:
		particle_pos
		particle_vel
		num_particles

	MAC Cells:
		cell_rows
		cell_cols
		cell_dim
		cell_types				# Type of each cell. {CellType::FLUID, CellType::AIR, CellType::SOLID}
		cell_velocities			# Vx (sampled at the left edge of the cell)
								# Vy (sampled at the bottom edge of the cell)

		prev_cell_velocities

	A 2d vector to store the normalizing terms for cell velicities when to_cell is true

	*/


	if (to_cell)
	{
		// TODO: Reset velocities and update cell types
		for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
		{
			cell_velocities[i][j] = Eigen::Vector2f::Zero();
			prev_cell_velocities[i][j] = Eigen::Vector2f::Zero();

			if (cell_types[i][j] != CellType::SOLID) {
				cell_types[i][j] = CellType::AIR;
			}
		}
		// assign cell types!!
		for (int i = 0; i < num_particles; ++i)
		{
			// Eigen::Vector2f& pos = particle_pos[i];
			int cell_x = static_cast<int>(std::floor(particle_pos[i].x() / cell_dim));
			int cell_y = static_cast<int>(std::floor(particle_pos[i].y() / cell_dim));

			if (cell_x >= 0 && cell_x < cell_cols && cell_y >= 0 && cell_y < cell_rows) {
				if (cell_types[cell_y][cell_x] != CellType::SOLID) {
					cell_types[cell_y][cell_x] = CellType::FLUID;
				}
			}
		}
	}

	
	// If you find iterating over x and y components unintuitive, feel free to change the structure of the code
	for (int component = 0; component < 2; ++component)
	{
		// Might need some setup here
		std::vector<std::vector<float>> weight_accumulator(cell_rows, std::vector<float>(cell_cols, 0.0f));

		for (int i = 0; i < num_particles; ++i)
		{
			// TODO: Bilinear interpolation on staggered grid
			Eigen::Vector2f& pos = particle_pos[i];
			Eigen::Vector2f& vel = particle_vel[i];

			float base_x = pos.x() / cell_dim;
			float base_y = pos.y() / cell_dim;
			component == 0 ? base_y-=0.5f : base_x-=0.5f;
			int cell_x = static_cast<int>(std::floor(base_x));
			int cell_y = static_cast<int>(std::floor(base_y));
			float dx = base_x - cell_x;
			float dy = base_y - cell_y;

			// weights
			float w1 = (1 - dx) * (1 - dy);
			float w2 = dx * (1 - dy);
			float w3 = dx * dy;
			float w4 = (1 - dx) * dy;

			if(cell_x < 0 || cell_x + 1 >= cell_cols || cell_y < 0 || cell_y + 1>= cell_rows) continue;

			if (to_cell)
			{
				// TODO: Transfer particle velocities to cells and store weights
				if (component == 0) { // x component
					
					cell_velocities[cell_y][cell_x].x() += w1 * vel.x();
					weight_accumulator[cell_y][cell_x] += w1;

					cell_velocities[cell_y][cell_x + 1].x() += w2 * vel.x();
					weight_accumulator[cell_y][cell_x + 1] += w2;
					
					cell_velocities[cell_y + 1][cell_x + 1].x() += w3 * vel.x();
					weight_accumulator[cell_y + 1][cell_x + 1] += w3;
					
					cell_velocities[cell_y + 1][cell_x].x() += w4 * vel.x();
					weight_accumulator[cell_y + 1][cell_x] += w4;
				}
				else { // y component
					cell_velocities[cell_y][cell_x].y() += w1 * vel.y();
					weight_accumulator[cell_y][cell_x] += w1;
					
					cell_velocities[cell_y][cell_x + 1].y() += w2 * vel.y();
					weight_accumulator[cell_y][cell_x + 1] += w2;
					
					cell_velocities[cell_y + 1][cell_x + 1].y() += w3 * vel.y();
					weight_accumulator[cell_y + 1][cell_x + 1] += w3;
					
					cell_velocities[cell_y + 1][cell_x].y() += w4 * vel.y();
					weight_accumulator[cell_y + 1][cell_x] += w4;
				}
			}
			else
			{
				// TODO: Transfer valid cell velocities back to the particles using a mixture of PIC and FLIP
				Eigen::Vector2f pic_velocity = Eigen::Vector2f::Zero();
				Eigen::Vector2f flip_velocity = Eigen::Vector2f::Zero();

				if (component == 0) { // x component
					float weight_sum_x = 0.0f;
					float vx1 = 0.0f, vx2 = 0.0f, vx3 = 0.0f, vx4 = 0.0f;
					float vx1_prev = 0.0f, vx2_prev = 0.0f, vx3_prev = 0.0f, vx4_prev = 0.0f;
					
					if ((cell_types[cell_y][cell_x] == CellType::FLUID) || 
						(cell_x > 0 && cell_types[cell_y][cell_x - 1] == CellType::FLUID)) {
						vx1 = cell_velocities[cell_y][cell_x].x();
						vx1_prev = prev_cell_velocities[cell_y][cell_x].x();	
						weight_sum_x += w1;
					}
					
					if ((cell_types[cell_y][cell_x + 1] == CellType::FLUID) || 
						(cell_types[cell_y][cell_x] == CellType::FLUID)) {
						vx2 = cell_velocities[cell_y][cell_x + 1].x();
						vx2_prev = prev_cell_velocities[cell_y][cell_x + 1].x();
						weight_sum_x += w2;
					}
					
					if ((cell_types[cell_y + 1][cell_x + 1] == CellType::FLUID) || 
						(cell_types[cell_y + 1][cell_x] == CellType::FLUID)) {
						vx3 = cell_velocities[cell_y + 1][cell_x + 1].x();
						vx3_prev = prev_cell_velocities[cell_y + 1][cell_x + 1].x();
						weight_sum_x += w3;
					}
					
					if ((cell_types[cell_y + 1][cell_x] == CellType::FLUID) || 
						(cell_x > 0 && cell_types[cell_y + 1][cell_x - 1] == CellType::FLUID)) {
						vx4 = cell_velocities[cell_y + 1][cell_x].x();
						vx4_prev = prev_cell_velocities[cell_y + 1][cell_x].x();
						weight_sum_x += w4;
					}
					
					// PIC
					if (weight_sum_x > 0.0f) {
						pic_velocity.x() = (w1 * vx1 + w2 * vx2 + w3 * vx3 + w4 * vx4) / weight_sum_x;
					}
					// FLIP
					 if (weight_sum_x > 0.0f) {
						flip_velocity.x() = vel.x() + 
							(w1 * (vx1 - vx1_prev) + w2 * (vx2 - vx2_prev) + w3 * (vx3 - vx3_prev) + w4 * (vx4 - vx4_prev)) / weight_sum_x;
					}
					// TOTAL
					vel.x() = (1 - flip_ratio) * pic_velocity.x() + flip_ratio * flip_velocity.x();
				}
				else { // y component
					float weight_sum_y = 0.0f;
					float vy1 = 0.0f, vy2 = 0.0f, vy3 = 0.0f, vy4 = 0.0f;
					float vy1_prev = 0.0f, vy2_prev = 0.0f, vy3_prev = 0.0f, vy4_prev = 0.0f;

					if ((cell_types[cell_y][cell_x] == CellType::FLUID) || 
						(cell_y > 0 && cell_types[cell_y - 1][cell_x] == CellType::FLUID)) {
						vy1 = cell_velocities[cell_y][cell_x].y();
						vy1_prev = prev_cell_velocities[cell_y][cell_x].y();
						weight_sum_y += w1;
					}
					
					if ((cell_types[cell_y][cell_x + 1] == CellType::FLUID) || 
						(cell_y > 0 && cell_types[cell_y - 1][cell_x + 1] == CellType::FLUID)) {
						vy2 = cell_velocities[cell_y][cell_x + 1].y();
						vy2_prev = prev_cell_velocities[cell_y][cell_x + 1].y();
						weight_sum_y += w2;
					}
					
					if ((cell_types[cell_y + 1][cell_x + 1] == CellType::FLUID) || 
						(cell_types[cell_y][cell_x + 1] == CellType::FLUID)) {
						vy3 = cell_velocities[cell_y + 1][cell_x + 1].y();
						vy3_prev = prev_cell_velocities[cell_y + 1][cell_x + 1].y();
						weight_sum_y += w3;
					}

					if ((cell_types[cell_y + 1][cell_x] == CellType::FLUID) || 
						(cell_types[cell_y][cell_x] == CellType::FLUID)) {
						vy4 = cell_velocities[cell_y + 1][cell_x].y();
						vy4_prev = prev_cell_velocities[cell_y + 1][cell_x].y();
						weight_sum_y += w4;
					}
				
					// PIC
					if (weight_sum_y > 0.0f) {
						pic_velocity.y() = (w1 * vy1 + w2 * vy2 + w3 * vy3 + w4 * vy4) / weight_sum_y;
					}
					// FLIP
					 if (weight_sum_y > 0.0f) {
						flip_velocity.y() = vel.y() +
							(w1 * (vy1 - vy1_prev) + w2 * (vy2 - vy2_prev) + w3 * (vy3 - vy3_prev) + w4 * (vy4 - vy4_prev)) / weight_sum_y;
					}
					// TOTAL
					vel.y() = (1 - flip_ratio) * pic_velocity.y() + flip_ratio * flip_velocity.y();
				}
			}
		}

		if (to_cell)
		{
			// TODO: Normalize cell velocities and store a backup in prev_velocities.
			for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
			{
				if (component == 0){
					if(weight_accumulator[i][j] > 0.0f){
						cell_velocities[i][j].x() /= weight_accumulator[i][j];
					}
					prev_cell_velocities[i][j].x() = cell_velocities[i][j].x();
				}
				else {
					if(weight_accumulator[i][j] > 0.0f){
						cell_velocities[i][j].y() /= weight_accumulator[i][j];
					}
					prev_cell_velocities[i][j].y() = cell_velocities[i][j].y();
				}
			}
		}
	}
}

void Fluid::updateDensity()
{
	/*TODO*/
	/*
	Here we will update the cell densities, which will be used to determine if a cell is overly compressed in the next step.
	The density of each cell is sampled at the center, and it uses the same concept as the transferVelocities function. 

	1.	Perform bilinear interpolation on each particle, each particle will contribute 1 * w_n to each of the sampling points.

	2.	Set the the resting density of water cells
		
		particle rest density = sum of all water cell densities / number of water cells

		This only needs to be done once, so simply test if particle_rest_density is 0 (the initial value).
		The particles are initialized to be tightly packed, so if you changed the layout of the particles, this step will have incorrect behavior.


	Some of the variables you might need, but not limited to:

	particle_rest_density		# updated only once at the first pass (test if particle_rest_density == 0)

	MAC cells:
		cell_densities
		cell_types
		cell_row
		cell_cols
		cell_dim

	Particle Information:
		num_particles
		particle_pos
	*/

	for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
		cell_densities[i][j] = 0.0f;

	
	for (int i = 0; i < num_particles; ++i)
	{
		// TODO: Perform bilinear interpolation
		Eigen::Vector2f& pos = particle_pos[i];
		float base_x = (pos.x() - 0.5f * cell_dim) / cell_dim;
		float base_y = (pos.y() - 0.5f * cell_dim) / cell_dim;
		int cell_x = static_cast<int>(std::floor(base_x));
		int cell_y = static_cast<int>(std::floor(base_y));
		float dx = base_x - cell_x;
		float dy = base_y - cell_y;

		// weights
		float w1 = (1 - dx) * (1 - dy);
		float w2 = dx * (1 - dy);
		float w3 = dx * dy;
		float w4 = (1 - dx) * dy;

		if (cell_x >= 0 && cell_x + 1 < cell_cols && cell_y >= 0 && cell_y + 1 < cell_rows) {
			cell_densities[cell_y][cell_x] += w1;
			cell_densities[cell_y][cell_x + 1] += w2;
			cell_densities[cell_y + 1][cell_x + 1] += w3;
			cell_densities[cell_y + 1][cell_x] += w4;
		}
	}

	if (particle_rest_density == 0.0) 
	{
		// TODO: Calculate resting particle densities in fluid cells.
		float total_density = 0.0f;
		int num_fluid_cells = 0;
		for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
		{
			if (cell_types[i][j] == CellType::FLUID) {
				total_density += cell_densities[i][j];
				num_fluid_cells++;
			}
		}
		if (num_fluid_cells > 0) {
			particle_rest_density = total_density / num_fluid_cells;
		}
	}

}

void Fluid::solveIncompressibility()
{
	/*TODO*/
	/*
	We will implement an iterative solver to enforce incompressibility on the cells.
	Here's what you will need to do:

	For each iteration, iterate over all fluid cells and perform the following:

		Sum up the number of potential flow directions (number of non-solid neighbor cells).

		Calculate the divergence of a cell using cell velocities.

		Calculate how much the cell density differs from the rest density and if the cell is overly compressed,
		we treat it as having more inflow by adding a penalty of stiffness_coefficient * compression to the divergence 
		to promote the expansion of fluids within the cell.

		The divergence needs to be corrected using the available flow directions, so the correction magnitude will be:
			divergence / normalizing_term * over_relaxation
		
		Apply the correction to cell velocities to try and solve for divergence = 0.

	Some of the variables you might need, but not limited to:

	particle_rest_density		
	density_correction			# Boolean that controls whether density correction is applied
	stiffness_coefficient		# Constant you can set with sliders
	over_relaxation				# Constant you can set with sliders

	MAC Cells:
		cell_rows
		cell_cols
		cell_types
		cell_velocities
		cell_densities
		
	*/

	// TODO: Store velocities
	// std::vector<std::vector<Eigen::Vector2f>> original_velocities = cell_velocities;

	for (int iter = 0; iter < iterations; ++iter)
	{
		for (int i = 1; i < cell_rows - 1; ++i) for (int j = 1; j < cell_cols - 1; ++j)
		{
			// TODO: Calculate divergence of fluid cells
			if (cell_types[i][j] != CellType::FLUID) continue;
			float divergence = 0.0f;
			if (cell_types[i][j - 1] != CellType::SOLID) divergence -= cell_velocities[i][j].x();
			if (cell_types[i][j + 1] != CellType::SOLID) divergence += cell_velocities[i][j+1].x();
			if (cell_types[i - 1][j] != CellType::SOLID) divergence -= cell_velocities[i][j].y();
			if (cell_types[i + 1][j] != CellType::SOLID) divergence += cell_velocities[i+1][j].y();
			// float divergence = cell_velocities[i][j+1].x() + cell_velocities[i+1][j].y() - cell_velocities[i][j].x() - cell_velocities[i][j].y();
			// TODO: Add bias to ouflow if density_correction is true
			if (density_correction) {
				float density_diff = cell_densities[i][j] - particle_rest_density;
				if (density_diff > 0.0f) {
					divergence -= stiffness_coefficient * density_diff;
				}
			}
			// TODO: Correct the cell velocities
			int num_directions = 0;
			if (cell_types[i - 1][j] != CellType::SOLID) num_directions++;
			if (cell_types[i + 1][j] != CellType::SOLID) num_directions++;
			if (cell_types[i][j - 1] != CellType::SOLID) num_directions++;
			if (cell_types[i][j + 1] != CellType::SOLID) num_directions++;
			if (num_directions > 0) {
				float correction = over_relaxation * divergence / num_directions;
				if (cell_types[i][j - 1] != CellType::SOLID) cell_velocities[i][j].x() += correction;

				if (cell_types[i][j + 1] != CellType::SOLID) cell_velocities[i][j + 1].x() -= correction;

				if (cell_types[i - 1][j] != CellType::SOLID) cell_velocities[i][j].y() += correction;

				if (cell_types[i + 1][j] != CellType::SOLID) cell_velocities[i + 1][j].y() -= correction;
			}
		}
	}
}

Fluid::Fluid(int num_particles, float radius, float obstacle_radius, float flip_ratio, float cell_dim, float relaxation_cell_dim, int iterations, float viewport_w, 
	float viewport_h, float dt, float gravity, float stiffness, bool density_correction, float over_relaxation)
	: num_particles(num_particles), particle_radius(radius), cell_dim(cell_dim), relaxation_cell_dim(relaxation_cell_dim), viewport_w(viewport_w), viewport_h(viewport_h), 
	dt(dt), gravity(gravity), obstacle_x(viewport_w), obstacle_y(viewport_h), obstacle_vx(0.0f), obstacle_vy(0.0f), obstacle_r(obstacle_radius), iterations(iterations), 
	flip_ratio(flip_ratio), stiffness_coefficient(stiffness), density_correction(density_correction), over_relaxation(over_relaxation)
{
	cell_rows = static_cast<int>(viewport_h / cell_dim);
	cell_cols = static_cast<int>(viewport_w / cell_dim);

	relaxation_cell_rows = static_cast<int>(viewport_h / relaxation_cell_dim);
	relaxation_cell_cols = static_cast<int>(viewport_w / relaxation_cell_dim);

	relaxation_cell_particle_ids = std::vector<std::vector<std::vector<int>>>(relaxation_cell_rows, std::vector<std::vector<int>>(relaxation_cell_cols));
	cell_types = std::vector<std::vector<CellType>>(cell_rows, std::vector<CellType>(cell_cols, CellType::AIR));
	cell_velocities = std::vector<std::vector<Eigen::Vector2f>>(cell_rows, std::vector<Eigen::Vector2f>(cell_cols, Eigen::Vector2f::Zero()));
	prev_cell_velocities = std::vector<std::vector<Eigen::Vector2f>>(cell_rows, std::vector<Eigen::Vector2f>(cell_cols, Eigen::Vector2f::Zero()));
	cell_colors = std::vector<Eigen::Vector3f>(cell_rows * cell_cols, Eigen::Vector3f::Zero());
	cell_densities = std::vector<std::vector<float>>(cell_rows, std::vector<float>(cell_cols, 0.0f));
	cell_centers_rendering.reserve(cell_rows * cell_cols);

	particle_pos = std::vector<Eigen::Vector2f>(num_particles, Eigen::Vector2f::Zero());
	particle_vel = std::vector<Eigen::Vector2f>(num_particles, Eigen::Vector2f::Zero());
	particle_colors = std::vector<Eigen::Vector3f>(num_particles, Eigen::Vector3f::Ones());

	for (int i = 1; i < cell_rows; ++i)
	{
		cell_types[i][0] = CellType::SOLID;
		cell_types[i][cell_cols - 1] = CellType::SOLID;
	}

	for (int j = 0; j < cell_cols; ++j)
	{
		cell_types[0][j] = CellType::SOLID;
		cell_types[cell_rows - 1][j] = CellType::SOLID;
	}

	// You can modify the way particle positions are initialized,
	// but the layout of the particles should be the same since it is used for rest density calculation.
	float x = cell_dim * 10;
	float y = cell_dim * 10;
	float x_max = 0.5f * viewport_w;

	float dx = 2 * radius;
	float dy = sqrt(3.0) / 2.0 * dx;

	bool stagger = false;

	for (int i = 0; i < num_particles; i++) 
	{
		particle_pos[i] = Eigen::Vector2f(x, y);
		x += dx;
		if (x > x_max)
		{
			y += dy;
			x = stagger ? cell_dim * 10 + radius : cell_dim * 10;
			stagger = !stagger;
		}
	}
}

void Fluid::update(int render_option, bool lmbDown)
{

	// physics
	integrate();
	particleRelaxation();
	handleCollisions(lmbDown);
	transferVelocities(true);
	updateDensity();
	solveIncompressibility();
	transferVelocities(false);

	
	// graphics
	switch (render_option)
	{
	case 0: // cells only
		updateCellColors();
		updateCellColorBuffers();
		renderCells();
		break;
	case 1: // both
		updateCellColors();
		updateCellColorBuffers();
		renderCells();
	case 2: // particles only
		updateParticleColors();
		updateParticleBuffers();
		renderParticles();
		break;
	default:
		return;
	}
	
	renderObstacle();
}

void Fluid::setObstacle(float mouse_x, float mouse_y, float mouse_vx, float mouse_vy)
{
	obstacle_x = mouse_x;
	obstacle_y = viewport_h - mouse_y;
	obstacle_vx = mouse_vx;
	obstacle_vy = -mouse_vy;
}

void Fluid::updateParticleColors()
{
	for (int i = 0; i < num_particles; ++i)
	{
		float vel = particle_vel[i].norm();
		float mult = std::min(vel, 50.0f) / 50.0f;
		particle_colors[i] = Eigen::Vector3f(mult, mult, 1.0f);
	}
}

void Fluid::updateCellColors()
{
	float max_density_estimate = cell_dim / particle_radius;
	for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
	{
		int ind = i * cell_cols + j;
		if (cell_types[i][j] == CellType::SOLID)
			cell_colors[ind] = Eigen::Vector3f(0.5f, 0.5f, 0.5f);
		else if (cell_types[i][j] == CellType::FLUID)
		{
			float mult = std::min(cell_densities[i][j], max_density_estimate) / max_density_estimate;
			
			cell_colors[ind] = Eigen::Vector3f(mult, mult, 1.0f);
		}
		else
			cell_colors[ind] = Eigen::Vector3f::Zero();
	}
}

void Fluid::setupRendering(GLuint pointShaders, GLuint obstacleShaders)
{
	// particles
	point_shaders = pointShaders;
	obstacle_shaders = obstacleShaders;

	glGenVertexArrays(1, &particles_vao);
	glBindVertexArray(particles_vao);

	glGenBuffers(1, &particles_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, particles_vbo);
	glBufferData(GL_ARRAY_BUFFER, num_particles * sizeof(Eigen::Vector2f),
		particle_pos.data(), GL_DYNAMIC_DRAW);

	GLint particleAttrPositionLoc = glGetAttribLocation(point_shaders, "attrPosition");
	if (particleAttrPositionLoc != -1) {
		glVertexAttribPointer(particleAttrPositionLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(particleAttrPositionLoc);
	}

	glGenBuffers(1, &particle_colors_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, particle_colors_vbo);
	glBufferData(GL_ARRAY_BUFFER, num_particles * sizeof(Eigen::Vector3f),
		particle_colors.data(), GL_DYNAMIC_DRAW);

	GLint particleAttrColorLoc = glGetAttribLocation(point_shaders, "attrColor");
	if (particleAttrColorLoc != -1) {
		glVertexAttribPointer(particleAttrColorLoc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glEnableVertexAttribArray(particleAttrColorLoc);
	}

	glBindVertexArray(0);

	// cells
	int num_cells = cell_rows * cell_cols;
	glGenVertexArrays(1, &cells_vao);
	glBindVertexArray(cells_vao);

	glGenBuffers(1, &cells_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, cells_vbo);

	for (int i = 0; i < cell_rows; ++i)
		for (int j = 0; j < cell_cols; ++j)
			cell_centers_rendering[i * cell_cols + j] = Eigen::Vector2f(j + 0.5, i + 0.5) * cell_dim;

	glBufferData(GL_ARRAY_BUFFER, num_cells * sizeof(Eigen::Vector2f),
		cell_centers_rendering.data(), GL_DYNAMIC_DRAW);

	GLint cellAttrPositionLoc = glGetAttribLocation(point_shaders, "attrPosition");

	glVertexAttribPointer(cellAttrPositionLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(cellAttrPositionLoc);

	glGenBuffers(1, &cell_colors_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, cell_colors_vbo);

	glBufferData(GL_ARRAY_BUFFER, num_cells * sizeof(Eigen::Vector3f),
		cell_colors.data(), GL_DYNAMIC_DRAW);

	GLint cellAttrColorLoc = glGetAttribLocation(point_shaders, "attrColor");

	glVertexAttribPointer(cellAttrColorLoc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(cellAttrColorLoc);

	glBindVertexArray(0);

	// obstacle dummy vao
	glGenVertexArrays(1, &obstacle_vao);
}

void Fluid::updateParticleBuffers()
{
	glBindBuffer(GL_ARRAY_BUFFER, particles_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, num_particles * sizeof(Eigen::Vector2f), particle_pos.data());

	glBindBuffer(GL_ARRAY_BUFFER, particle_colors_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, num_particles * sizeof(Eigen::Vector3f), particle_colors.data());
}

void Fluid::renderParticles()
{
	glUseProgram(point_shaders);

	glUniform2f(glGetUniformLocation(point_shaders, "domainSize"), viewport_w, viewport_h);
	glUniform1f(glGetUniformLocation(point_shaders, "pointSize"), 2.0f * particle_radius);
	glUniform1f(glGetUniformLocation(point_shaders, "drawDisk"), 1.0f);

	glBindVertexArray(particles_vao);
	glDrawArrays(GL_POINTS, 0, num_particles);

	glBindVertexArray(0);
}

void Fluid::updateCellColorBuffers()
{
	glBindBuffer(GL_ARRAY_BUFFER, cell_colors_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, cell_rows * cell_cols * sizeof(Eigen::Vector3f), cell_colors.data());

}
void Fluid::renderCells()
{
	glUseProgram(point_shaders);

	glUniform2f(glGetUniformLocation(point_shaders, "domainSize"), viewport_w, viewport_h);
	glUniform1f(glGetUniformLocation(point_shaders, "pointSize"), cell_dim);
	glUniform1f(glGetUniformLocation(point_shaders, "drawDisk"), 0.0f);

	glBindVertexArray(cells_vao);
	glDrawArrays(GL_POINTS, 0, cell_rows * cell_cols);

	glBindVertexArray(0);
}

void Fluid::renderObstacle()
{
	glUseProgram(obstacle_shaders);

	glBindVertexArray(obstacle_vao);

	glUniform2f(glGetUniformLocation(obstacle_shaders, "domainSize"), viewport_w, viewport_h);
	glUniform2f(glGetUniformLocation(obstacle_shaders, "attrPosition"), obstacle_x, obstacle_y);
	glUniform1f(glGetUniformLocation(obstacle_shaders, "pointSize"), 2 * obstacle_r);
	
	glDrawArrays(GL_POINTS, 0, 1);
	glBindVertexArray(0);
}



