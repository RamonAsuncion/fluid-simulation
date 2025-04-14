[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/0Ct0PJdd)
# Quest 4: The Enchanted Symphony of Motion

## Objectives
As a budding wizard of the Computer Graphics Wizard Academy, your next quest is to summon and control dynamic visual phenomena. In this quest, you will choose between conjuring a particle system or crafting a mass-spring system. Both require you to weave powerful spells in WebGPU, combining knowledge of physics and animation to create interactive and enchanting visual effects.

## Quest Overview
The Enchanted Symphony of Motion requires apprentices to: 

1. Summon Particles or Springs: use WebGPU to create mesmerizing visual effects like fire, smoke, sparkles, or elastic behaviors such as cloth waving or pendulum oscillations.
2. Interactive Wizardry: implement keyboard and mouse spells to influence the simulation in real time.
3. Arcane Performance: utilize GPU resources effectively to ensure your visual spells perform smoothly at an interactive frame rate.
   
## Feature Points
Earn feature point listed below is worth **1 point** except the race condition one. The total possible points for this quest are **10+** for the entire quest. You only need 60 points from 10 quests. Do as many as you want. On average, you should get 6 points per quest.

**IMPORTANT**: In order to receive points, you **must** implement a different visual effect from the one you have learned in the scrolls, or show that you significantly improved the particle/mass-spring system.

### Particle System
- **1 point**: Setup a particle system physics in compute shaders. (i.e. initial particles, emitters etc.)
- **1 point**: Introduce nature forces such as gravity, wind, central attraction force, to influence particle motion.
- **1 point**: Add lifespan for particles (i.e. removing them as they age) and keep a max bound of the number of particles.
- **1 point**: Update particle positions based on velocity and acceleration (forces).
- **1 point**: Implement circular boundary condition. Particle goes off screen will wrap around to the other side.
- **1 point**: Use mouse/keyboard interaction to spawn or attract particles.
- **1 point**: Map particles with textures to create visual apperances.
- **1 point**: Create a mesmerizing effect such as fireworks, smoke, or water etc.

### Mass-Spring System
- **1 point**: Setup a mass-spring system physics in compute shaders.
- **1 point**: Introduce nature forces such as gravity, wind to influence the tension/compression between particles.
- **1 point**: Add spring force and damping to stabilaize the spring oscillations.
- **1 point**: Update particle positions based on velocity and acceleration (spring forces and damping).
- **1 point**: Pin/fix specific particles to act as anchors, preventing the system from collapsing.
- **1 point**: Use mouse/keyboard interaction to select particles and apply external forces to them.
- **1 point**: Map surfaces with textures to create visual apperances.
- **1 point**: Simulate a realistic cloth, jelly-like stwructure, etc.
- **2 points**: Implement a scheme that avoids **race conditions**. e.g. graph coloring, or coverting floats to ints and use atomic operations.

### Common
- **1 point**: Apply non-linear interpolation.
- **1 point**: Compose a magical scene that has at least two different effects. You do not need to implement both systems. You can reuse the same system to design a secondary effect.
- **1 point**: Run at real time even with a large amount of particles (e.g. $\ge$ 10,000).
- **1 point**: Win the weekly quest contest. On Monday, in class, you will vote for the best craft of the week. **If you win, you earn an extra point!**
- **? points**: Anything that impresses me (stunning visual effects using the particle system such as water fall, hair simulation using the mass-spring system, optimal graph positioning using the mass-spring system etc.). **Write them clearly in your description to convince me they are impressive!**

### Submission Instructions
- Fork this quest repository to begin your quest and ensure you commit your code there before the due date.
- Ensure your Arcane Portal is public and contains all relevant and **obfuscated** code. I recommend you update the Arcane Portal every time you complete a quest.
- Your Arcane Portal site should display the visual effect(s) and include any required descriptions.

### Resources
- Refer to **[Scroll 9](https://eg.bucknell.edu/~scl019/Courses/CGSP25/scroll9.php)** for 2D particle systems.
- Consult **[Scroll 10](https://eg.bucknell.edu/~scl019/Courses/CGSP25/scroll10.php)** for more particle system eamples
- Check **[Scroll 11](https://eg.bucknell.edu/~scl019/Courses/CGSP25/scroll11.php)** for 2D mass-spring systems.
- [WGSL documentation](https://www.w3.org/TR/WGSL/).
- Reach out and ask questions on [Piazza](https://piazza.com/bucknell/spring2025/csci379) if you are stuck!
