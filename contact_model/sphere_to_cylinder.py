from utils.vector_algebra import unit_vector, dot_product
from utils.geometry import find_point_projection_along_a_line, get_distance_between_edges

# type hint imports
from numpy import ndarray
from opensim import Model, State, Body, Frame, Marker, Vec3, Ground

def compute_x_and_x_dot(model: Model, sphere_loc: ndarray, cylinder_vel: ndarray, clavicle: Body, sphere_r: float,
                        cylinder_r: float, s: State):

    """
    This function aims at computing the penetration between a sphere and a cylinder and its 1st time derivative.

    :param model: OpenSim Model object.
    :param sphere_loc: location of the sphere origin in its parent Body reference frame.
    :param cylinder_vel: cylinder body velocity in the global reference frame.
    :param clavicle: clavicle body.
    :param sphere_r: radius of the sphere.
    :param cylinder_r: radius of the cylinder.
    :param s: OpenSim state
    :return:
    """

    # ground
    ground: Ground = model.getGround()

    # convert sphere location to Vec3 object
    sphere_loc_vec3: Vec3 = Vec3.createFromMat(sphere_loc.squeeze(-1))

    # get clavicle Base Frame
    clavicleFrame: Frame = clavicle.findBaseFrame()

    # find sphere CoM in global reference frame
    sphere_com: Vec3 = clavicleFrame.findStationLocationInGround(s, sphere_loc_vec3)

    # get cylinder frame
    cylinderBody: Body = model.getBodySet().get("punching_bag")
    cylinderFrame: Frame = cylinderBody.findBaseFrame()

    # ------------------------------------------------------------------------------------------------------------------
    # compute the position of the centre of the top and bottom surfaces of the cylinder
    # ------------------------------------------------------------------------------------------------------------------
    # find top centre
    cylinder_top_marker: Marker = model.getMarkerSet().get("cylinder_top")
    cylinder_top: Vec3 = cylinder_top_marker.findLocationInFrame(s, model.getGround())

    # find bottom centre
    cylinder_bottom_marker: Marker = model.getMarkerSet().get("cylinder_bottom")
    cylinder_bottom: Vec3 = cylinder_bottom_marker.findLocationInFrame(s, model.getGround())


    # find projection of the sphere CoM onto the cylinder longitudinal axis --> this will be a relative position vector,
    # expressed with respect to the bottom surface center of the cylinder.
    Q_relative: ndarray = find_point_projection_along_a_line(sphere_com.to_numpy()[:, None],
                                                          cylinder_bottom.to_numpy()[:, None],
                                                          cylinder_top.to_numpy()[:, None])

    # Now we add the cylinder bottom center vector to express Q in the Global reference frame.
    Q_global: ndarray = Q_relative + cylinder_bottom.to_numpy()[:, None]

    # compute directional vector of cylinder velocity
    cylinder_motion_direction: ndarray = unit_vector(cylinder_vel)


    # key function: here we compute the indentation between the 2 surfaces. This quantity is a scalar that consists of
    # the 3D relative position vector between the 2 surface edges, projected onto the directional vector that points from
    # the center of the sphere to Q. This directional vector is outputted as "n".
    d, cylinder_edge, sphere_edge, n = get_distance_between_edges(Q_global, sphere_com.to_numpy()[:, None],
                                                               sphere_r, cylinder_r, cylinder_motion_direction)

    # ------------------------------------------------------------------------------------------------------------------
    # Compute rate of change over time of indentation as the relative velocity of the 2 edge points
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

    # local vectors: express the point on the edges in their respective reference frame.
    cylinderEdge_station: Vec3 = ground.findStationLocationInAnotherFrame(s, cylinderEdge_loc, cylinderFrame)
    sphereEdge_station: Vec3 = ground.findStationLocationInAnotherFrame(s, sphereEdge_loc, clavicleFrame)

    # find velocity of these 2 points in ground
    cylinderEdge_vel: Vec3 = cylinderFrame.findStationVelocityInGround(s, cylinderEdge_station)
    spheredEdge_vel: Vec3 = clavicleFrame.findStationVelocityInGround(s, sphereEdge_station)

    # compute relative vel
    vel: ndarray = cylinderEdge_vel.to_numpy() - spheredEdge_vel.to_numpy()

    # project vel onto normal unit vector
    vel_scalar: float = float(dot_product(vel, n))

    return d, vel_scalar
