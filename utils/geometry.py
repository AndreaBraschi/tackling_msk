from utils.vector_algebra import *
from typing import Tuple

# type hint imports
from numpy import ndarray


def find_point_projection_along_a_line(P0: ndarray, P1: ndarray, P2: ndarray) -> ndarray:
    """
        This function computes the projection of a 3D point onto a line.

        Args:
            P0: three-dimensional point that we want to project on the line.
            P1: 3D point that represent the bottom endpoint of the line.
            P2: 3D point that represent the top endpoint of the line.

        """

    # 1) find vector that goes from one point along the line to the point
    r: ndarray = find_vector_between_two_points(P1, P2)

    # 2) find direction of the vector
    r_hat: ndarray = unit_vector(r)

    # 3) compute vector that goes from one of the 2 points to the one that doesn't lie on the line
    q: ndarray = find_vector_between_two_points(P1, P0)

    # 4) project q onto r_hat to find the scalar value t that dictates how much r_hat will have to be sclaed uo or
    #    down to go from P1 to the poit projection of P0 onto the r-line
    t: ndarray = dot_product(q, r)

    # 5) compute the cosine of the angle between the 2 vectors
    cos_theta = t / (norm(q) * norm(r))

    # 6) project q onto r
    projection = q * cos_theta

    # 5) scale the unit vector to get the 3D coordinate of P0 projection onto r
    Q: ndarray = projection * r_hat

    if Q.shape[-2] < Q.shape[-1]:
        Q = np.transpose(Q)

    return Q

def get_distance_between_edges(Q: ndarray, sphere_com: ndarray, radius_sphere: float, radius_cylinder: float,
                               motion_direction: ndarray) -> Tuple[ndarray, ...]:

    """
    This function computes the distance between the edges of the sphere and the cylinder.

    :param Q: the point along the cylinder longitudinal axis that resulted to be the closest to the sphere.
    :param sphere_com: sphere centre of mass. 3 x 1 vector
    :param radius_sphere: sphere radius
    :param radius_cylinder: cylinder radius
    :param motion_direction: direction of motion of the cylinder
    :return:
    """


    # Here, for both geometries, we need to make sure that we go to the edge that corresponds to the direction where the
    # geometries are actually moving

    # distance vector (3 x 1) from Q to the spere centre of mass
    r: ndarray = sphere_com - Q
    # directional vector that points from Q at sphere_com
    r_hat: ndarray = unit_vector(r)
    # directional vector that points from sphere_com at Q
    r_hat_: ndarray = unit_vector(r_hat * (-1))

    print(f"motion direction: {motion_direction[0]}")
    print(f"r_hat: {r_hat[0]}")

    # here we have handle the 2 following cases
    # 1) when the sphere centre of mass hasn't gone over the cylinder edge: in this case we have the cylinder horizontal
    #    motion direction vector that is pointing in the same direction as the unit vector that points towards the sphere
    #    CoM
    if motion_direction[0] * r_hat[0] > 0:

        if motion_direction[0] > 0:
            cylinder_edge: ndarray = Q + (radius_cylinder * r_hat * (-1))
            sphere_edge: ndarray = sphere_com + (radius_sphere * r_hat_ * (-1))

        else:
            cylinder_edge: ndarray = Q + (radius_cylinder * r_hat)
            sphere_edge: ndarray = sphere_com + (radius_sphere * r_hat_)


    # 2) when the sphere centre of mass has gone over the cylinder edge: in this case we have the cylinder horizontal
    #    motion direction vector that is pointing in the opposite direction as the unit vector that points towards the
    #    sphere CoM.
    else:

        if motion_direction[0] > 0:
            cylinder_edge: ndarray = Q + (radius_cylinder * r_hat)
            sphere_edge: ndarray = sphere_com + (radius_sphere * r_hat_)

        else:
            cylinder_edge: ndarray = Q + (radius_cylinder * r_hat * (-1))
            sphere_edge: ndarray = sphere_com + (radius_sphere * r_hat_ * (-1))


    # distance vector (3 x 1) between the 2 edges in the global reference frame.
    d: ndarray = cylinder_edge - sphere_edge

    return d, cylinder_edge, sphere_edge