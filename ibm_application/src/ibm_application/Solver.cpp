#include "ibm_application/Solver.h"

void Solver::PerformStep(int steps)
{
	for (size_t i = 0; i < steps; i++)
	{
		if (m_time+m_dt > m_end_time)
		{
			break;
		}

		m_selected_scheme->Update(m_dt, m_cfl);
		m_time += m_dt;
	}
}