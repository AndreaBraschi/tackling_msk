from numpy import array
from utils.vector_algebra import unit_vector, norm
from utils.geometry import find_point_projection_along_a_line, get_distance_between_edges

# type hint imports
from numpy import ndarray
from opensim import State, Body, Frame, Rotation, Matrix

def compute_penetration(sphere_loc: ndarray, cylinder: Body, cylinder_vel: ndarray, clavicle: Body, s: State):

    """
    This function aims at computing the penetration between a sphere and a cylinder

    :param sphere_loc: location of the sphere origin in its parent Body reference frame.
    :param cylinder: cylinder body.
    :param cylinder_vel: cylinder body velocity in the global reference frame.
    :param clavicle: clavicle body.
    :param s: OpenSim state
    :return:
    """

    # cylinder half height
    half_h: float = 0.27
    # create a 3x1 vector that would represent the orientation of the vector that goes from the cylinder body CoM to
    # the centre of the top surface, when the cylinder is standing perfectly aligned with the OpenSim Y-axis.
    half_h_vec: ndarray = array([[0.0], [half_h], [0.0]])   # shape: [3, 1]

    # cylinder radius
    cylinder_r: float = 0.1702085

    # sphere radius
    sphere_r: float = 0.025

    # ------------------------------------------------------------------------------------------------------------------
    # compute pose of clavicle and cylinder
    # ------------------------------------------------------------------------------------------------------------------

    # ----- Get clavicle pose in global coordinate frame ----- #
    # get clavicle Base Frame
    clavicleFrame: Frame = clavicle.findBaseFrame()

    # get clavicle CoM in Global reference frame
    O_clavicle: ndarray = clavicleFrame.getPositionInGround(s).to_numpy()[:, None]  # 3 x 1 vector

    # get clavicle orientation in Global reference frame
    R_osim: Rotation = clavicleFrame.getRotationInGround(s)
    R_clavicle: ndarray = array([[R_osim.get(0, 0), R_osim.get(0, 1), R_osim.get(0, 2)],
                                [R_osim.get(1, 0), R_osim.get(1, 1), R_osim.get(1, 2)],
                                [R_osim.get(2, 0), R_osim.get(2, 1), R_osim.get(2, 2)]])


    # ----- Get cylinder pose in global coordinate frame ----- #
    # get cylinder Base Frame
    cylinderFrame: Frame = cylinder.findBaseFrame()

    # get cylinder CoM in Global reference frame
    O_cylinder: ndarray = cylinderFrame.getPositionInGround(s).to_numpy()[:, None]  # 3 x 1 vector

    # get cylinder orientation in Global reference frame
    R_osim: Rotation = cylinderFrame.getRotationInGround(s)
    R_cylinder: ndarray = array([[R_osim.get(0, 0), R_osim.get(0, 1), R_osim.get(0, 2)],
                                 [R_osim.get(1, 0), R_osim.get(1, 1), R_osim.get(1, 2)],
                                 [R_osim.get(2, 0), R_osim.get(2, 1), R_osim.get(2, 2)]])


    # ------------------------------------------------------------------------------------------------------------------
    # compute location of sphere geometric centre in the global reference frame
    # ------------------------------------------------------------------------------------------------------------------
    O_sphere: ndarray = R_clavicle @ sphere_loc + O_clavicle

    # ------------------------------------------------------------------------------------------------------------------
    # compute the position of the centre of the top and bottom surfaces of the cylinder
    # ------------------------------------------------------------------------------------------------------------------
    # rotate half height vector according to cylinder orientation
    oriented_half_h: ndarray = R_cylinder @ half_h_vec  # 3 x 1 vector
    # find top centre by adding the half height to the cylinder CoM
    cylinder_top: ndarray = O_cylinder + oriented_half_h   # 3 x 1 vector
    # find bottom centre by subtracting the half height from the cylinder CoM
    cylinder_bottom: ndarray = O_cylinder - oriented_half_h   # 3 x 1 vector

    # ------------------------------------------------------------------------------------------------------------------
    # compute the position of the centre of the top and bottom surfaces of the cylinder
    # ------------------------------------------------------------------------------------------------------------------
    # find projection of the sphere CoM onto the cylinder longitudinal axis
    Q_local: ndarray = find_point_projection_along_a_line(O_sphere, cylinder_bottom, cylinder_top)
    # we need it in the global reference frame
    Q_global: ndarray = Q_local + cylinder_bottom

    # compute directional vector of cylinder velocity
    cylinder_motion_direction: ndarray = unit_vector(cylinder_vel)

    d, cylinder_edge, sphere_edge = get_distance_between_edges(Q_global, O_sphere, sphere_r, cylinder_r, cylinder_motion_direction)
    if d[0] < 0:
        print(f"Sphere 1: contact with TOP cylinder detected")
        return norm(d)

    else:
        return None
