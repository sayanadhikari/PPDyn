Physics Model
=============
Molecular Dynamics simulation (MD) is a powerful method to study atomic and molecular processes.
In plasma physics, it has vast applications starting from fusion plasma to dusty plasma.
It solves the equation of motion using `Newton equations of motion <https://en.wikipedia.org/wiki/Equations_of_motion>`_ (or `Langevin dynamics <https://en.wikipedia.org/wiki/Langevin_dynamics>`_ which includes 
energy dissipation through an additional friction term):

.. math::

  \frac{\partial^{2}}{\partial {t}^2}\vec{r}_i = \frac{1}{m_i}\vec{f}_i
  \vec{f}_i=-\frac{\partial}{\partial r} V(\vec{r}_1(t), \vec{r}_i(t),~...,~\vec{r}_N(t))

where :math:`\vec{r}_i(t)` is the position of atom :math:`i` at time :math:`t` with mass :math:`m_i`, and :math:`V` is the interaction potential between all
:math:`N` involved species.
