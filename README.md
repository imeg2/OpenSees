# A family of minimum residual displacement methods as nonlinear solution schemes for equilibrium path following in structural mechanics
# Mostafa Salehi Ahmad-Abad | Ali Maghami | Morteza Ghalishooyan | Ahmad Shooshtari

# Abstract

In the realm of nonlinear structural mechanics, tracking load-displacement paths becomes intricate when approaching the limits of instability in both load and displacement. These complex scenarios often perplex conventional methodologies,posing challenges in accurately characterizing structural behavior. Over the past decades, multiple initiatives have addressed this intricacy, dating back to the late 1960s. Despite the integration of various algorithms into commercial finite element software, it is essential to acknowledge that no single algorithm universally solves all nonlinear structural problems. Moreover, many existing methods are sensitive to initial parameters like step-length and may struggle with large values. Drawing inspiration from Chan’s concepts of the late 1980s, this study introduces a family of minimum residual displacement methods for tracking equilibrium paths. Our formulation can be seen as a generalized framework that encompasses the existing method as a specific instance. Building on this foundation, we develop and apply straightforward techniques to control residual displacement in nodes, elements, or specific displacement components. This versatile family of methods offers diverse applications. We present a comprehensive library of methods implemented in OpenSees, enabling a comparison between our approach and the conventional minimum residual displacement method, as well as four well-established techniques (cylindrical arclength, generalized displacement, modified generalized displacement, and updated normal plane). We apply these methods to solve five intricate geometrically nonlinear problems involving truss, frame, and shell structures (both isotropic and orthotropic). The results highlight the efficacy and efficiency of our methodologies, emphasizing their advantages in two key domains: their enhanced ability to converge even in
highly complex behavioral scenarios and their capacity to reduce the required iteration count.

# OpenSees Source Code Repository

This git repository contains all revisions to OpenSees source code since Version 2.3.2.

Older revisions of the code are available upon request.

If you plan on collaborating or even using OpenSees as your base code it is highly recommended that
you FORK this repo to your own account and work on it there. We will not allow anybody to write to
this repo. Only PULL requests will be considered. To fork the repo click on the FORK at the top of this page.

For a brief outline on forking we suggest:
https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow

For a brief introduction to using your new FORK we suggest:
https://www.atlassian.com/git/tutorials/saving-changes

## Documentation
The documentation for OpenSees is being moved to a parellel github repo:
https://github.com/OpenSees/OpenSeesDocumentation

The documentation (in its present form) can be viewed in the browser using the following url:
https://OpenSees.github.io/OpenSeesDocumentation

## Build Instructions
Steps to build OpenSees on Windows, Linux, and Mac:
https://opensees.github.io/OpenSeesDocumentation/developer/build.html

## Modeling Questions
Issues related to modeling questions will be closed. Instead, post your modeling questions on the OpenSees 
message board or in the OpenSees Facebook group.
+ https://opensees.berkeley.edu/community
+ https://facebook.com/groups/opensees
