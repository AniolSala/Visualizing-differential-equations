# Visualising differential equations: Particle moving over the surface of a cone
---

## Obtaining the quations
To obtain the equations of a particle moving over a surface of a cone I used
[lagrange's equations](https://en.wikipedia.org/wiki/Lagrangian_mechanics). The coordinates that define the trajectory will be
$\varphi$ and $z$.

![3D coordinates](3d_coords.png)
![3D coordinates](cone_fig.png)

Why? As $\theta$ is constant on the sufrace of the cone, $tan(\theta) = \frac{d}{z}$ is also constant. 
Then we only need two coordinates to describe the motion of the particle (as expected: a surface is a two-dimensional object!). 
The equations to solve are the following:

\begin{equation}
\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{q_{i}}} - \frac{\partial \mathcal{L}}{\partial q_{i}} = 0 
\end{equation}

where $\mathcal{L} = T - V$ and $q_{1, 2} = z, \theta$. The potential energy is just $V = - m\;g\;z$. For the kinetic energy $T$
we have to write $(x,\;y\;z)$ in terms of $z$ and $\theta$:

\begin{equation}
\begin{split}
x =& \; d\;cos(\varphi) = tan(\theta)\;z\;cos(\varphi) \\
y =& \; d\;sin(\varphi) = tan(\theta)\;z\;sin(\varphi)\\
z =& \; z
\end{split}
\end{equation}

and their derivates:

\begin{equation}
\begin{split}
\dot{x} =& \; tan(\theta) \left( \dot{z} \; \cos(\varphi) - z\; sin(\varphi) \; \dot{\varphi} \right) \\ 
\dot{y} =& \; tan(\theta) \left( \dot{z} \; \sin(\varphi) + z\; cos(\varphi) \; \dot{\varphi} \right)\\
\dot{z} =& \; \dot{z}
\end{split}
\end{equation}

And that's all we need! The kinetic energy can be written as:

\begin{equation}
T = \frac{1}{2} m \; \left( \dot{x}^2 + \dot{y}^2 + \dot{z}^2 \right) =  \left(1 + tan^2(\theta)\right)\dot{z}^2 + 
tan^2(\theta) \, z^2\dot{\varphi}^2 
\end{equation}

Now, using the Lagrange's equations, we obtain (finally!) the differential equations of $z$ and $\theta$:

\begin{equation}
\begin{split}
\ddot{z} =& \frac{tan^2(\theta)\;z\,\dot{\theta}^2 - g}{1 + tan^2(\theta)} \\
\ddot{\theta} =& - \frac{2\, \dot{z}\, \dot{\theta}}{z}
\end{split}
\end{equation}

## Visualizing the equations

At least! At this point the only thing remaining is solving the equations and plot them. The way I have done this is using the
[Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Derivation_of_the_Runge%E2%80%93Kutta_fourth-order_method) of fourth order (see the code [here]()).

Here there are some animations I have done with Matplotlib, each one with different initial conditions.

![cone_image](particle_cone_1.gif)
![cone_image](particle_cone_2.gif)
![cone_image](particle_cone_3.gif)