from numpy import array
from utils.vector_algebra import unit_vector, norm
from utils.geometry import find_point_projection_along_a_line, get_distance_between_edges

# type hint imports
from numpy import ndarray
from opensim import Model, State, Body, Frame, Marker, Vec3, Ground

def compute_penetration(model: Model, sphere_loc: ndarray, cylinder_vel: ndarray, clavicle: Body, s: State):

    """
    This function aims at computing the penetration between a sphere and a cylinder and its 1st time derivative.

    :param sphere_loc: location of the sphere origin in its parent Body reference frame.
    :param cylinder_vel: cylinder body velocity in the global reference frame.
    :param clavicle: clavicle body.
    :param s: OpenSim state
    :return:
    """

    # ground
    ground: Ground = model.getGround()

    # cylinder half height
    half_h: float = 0.5386

    # cylinder radius
    cylinder_r: float = 0.1702085

    # sphere radius
    sphere_r: float = 0.025

    # convert sphere location to Vec3 object
    sphere_loc_vec3: Vec3 = Vec3.createFromMat(sphere_loc.squeeze(-1))

    # ------------------------------------------------------------------------------------------------------------------
    # compute pose of cylinder
    # ------------------------------------------------------------------------------------------------------------------
    # ----- Get cylinder pose in global coordinate frame ----- #
    # get clavicle Base Frame --> we need this to compute the cylinder top and bottom positions w.r.t the clavicle frame
    # using the OpenSim API.
    clavicleFrame: Frame = clavicle.findBaseFrame()

    # ------------------------------------------------------------------------------------------------------------------
    # compute the position of the centre of the top and bottom surfaces of the cylinder
    # ------------------------------------------------------------------------------------------------------------------
    # find top centre
    cylinder_top_marker: Marker = model.getMarkerSet().get("cylinder_top")
    cylinder_top: Vec3 = cylinder_top_marker.findLocationInFrame(s, model.getGround())

    # find bottom centre
    cylinder_bottom_marker: Marker = model.getMarkerSet().get("cylinder_bottom")
    cylinder_bottom: Vec3 = cylinder_bottom_marker.findLocationInFrame(s, model.getGround())

    # find sphere CoM in global reference frame
    sphere_com: Vec3 = clavicleFrame.findStationLocationInGround(s, sphere_loc_vec3)

    # ------------------------------------------------------------------------------------------------------------------
    # compute the position of the centre of the top and bottom surfaces of the cylinder
    # ------------------------------------------------------------------------------------------------------------------
    # find projection of the sphere CoM onto the cylinder longitudinal axis --> we're converting the Vec3 objects to
    # NumPy arrays.
    Q_local: ndarray = find_point_projection_along_a_line(sphere_com.to_numpy()[:, None], cylinder_bottom.to_numpy()[:, None], cylinder_top.to_numpy()[:, None])

    # Q_local is a relative positional vector that goes from the bottom of the cylinder to the projected point. However,
    # we need it in the reference frame of the clavicle.
    Q_clavicleFrame: ndarray = Q_local + cylinder_bottom.to_numpy()[:, None]

    # compute directional vector of cylinder velocity
    cylinder_motion_direction: ndarray = unit_vector(cylinder_vel)

    d, cylinder_edge, sphere_edge = get_distance_between_edges(Q_clavicleFrame, sphere_com.to_numpy()[:, None],
                                                               sphere_r, cylinder_r, cylinder_motion_direction)

    # ------------------------------------------------------------------------------------------------------------------
    # compute relative velocity of the 2 edge points
    # ------------------------------------------------------------------------------------------------------------------
    # create Vec3 objects for the 2 edges.
    # Vec3 accepts only 1D arrays, meaning that we need to check whether the numpy arrays are shape: (3, 1) --> 2D. If
    # they're, we need to squeeze the 2nd dimension.

    if cylinder_edge.ndim > 1:
        cylinder_edge = cylinder_edge.squeeze(-1)

    if sphere_edge.ndim > 1:
        sphere_edge = sphere_edge.squeeze(-1)

    # global vector
    cylinderEdge_loc: Vec3 = Vec3.createFromMat(cylinder_edge)
    sphereEdge_loc: Vec3 = Vec3.createFromMat(sphere_edge)

    # local vector
    cylinderEdge_station: Vec3 = ground.findStationLocationInAnotherFrame(s, cylinderEdge_loc, clavicleFrame)
    sphereEdge_station: Vec3 = ground.findStationLocationInAnotherFrame(s, sphereEdge_loc, clavicleFrame)

    # find velocity of these 2 points, which are expressed in the clavicle reference frame, in ground
    cylinderEdge_vel: Vec3 = clavicleFrame.findStationVelocityInGround(s, cylinderEdge_station)
    spheredEdge_vel: Vec3 = clavicleFrame.findStationVelocityInGround(s, sphereEdge_station)

    # compute relative vel
    vel: ndarray = cylinderEdge_vel.to_numpy() - spheredEdge_vel.to_numpy()

    return d, vel
