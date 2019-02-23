//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/obb.h"
#include "utilities/math_utils.h"
#include "processes/find_intersected_geometrical_objects_with_obb_process.h"

namespace Kratos
{

template<class TEntity>
FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::FindIntersectedGeometricalObjectsWithOBBProcess(
    ModelPart& rPart1,
    ModelPart& rPart2,
    const double BoundingBoxFactor
    ) : BaseType(rPart1, rPart2),
        mBoundingBoxFactor(BoundingBoxFactor)
{
    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::FindIntersectedGeometricalObjectsWithOBBProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : BaseType(rModel.GetModelPart(ThisParameters["first_model_part_name"].GetString()),
        rModel.GetModelPart(ThisParameters["second_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_first_model_part_name = ThisParameters["first_model_part_name"].GetString();
    const std::string& r_second_model_part_name = ThisParameters["second_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_first_model_part_name == "") << "first_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_second_model_part_name == "") << "second_model_part_name must be defined on parameters" << std::endl;

    // Setting the bounding box factor
    mBoundingBoxFactor = ThisParameters["bounding_box_factor"].GetDouble();

    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::HasIntersection(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    const std::size_t work_dim = rFirstGeometry.WorkingSpaceDimension(); // TODO: DOMAIN_SIZE should be considered for consistency with other implementations
    if (this->Is(BOUNDARY)) {
        if (work_dim == 2) {
            return this->HasIntersection2D(rFirstGeometry, rSecondGeometry);
        } else {
            return this->HasIntersection3D(rFirstGeometry, rSecondGeometry);
        }
    } else {
        if (work_dim == 2) {
            return BaseType::HasIntersection2D(rFirstGeometry, rSecondGeometry);
        } else {
            return BaseType::HasIntersection3D(rFirstGeometry, rSecondGeometry);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // The local dimensions
    const std::size_t local_dimension_1 = rFirstGeometry.LocalSpaceDimension();
    const std::size_t local_dimension_2 = rSecondGeometry.LocalSpaceDimension();

    // The edges
    PointerVector<GeometryType> r_edges_1 = rFirstGeometry.Edges();
    const std::size_t number_of_edges_1 = (local_dimension_1 < 2) ? 1 : r_edges_1.size();
    PointerVector<GeometryType> r_edges_2 = rSecondGeometry.Edges();
    const std::size_t number_of_edges_2 = (local_dimension_2 < 2) ? 1 : r_edges_2.size();

    // First geometry
    for (std::size_t i_1 = 0; i_1 < number_of_edges_1; ++i_1) {
        auto& r_edge_1 = *(r_edges_1.begin() + i_1);

        // Check the intersection of each edge of the object bounding box against the intersecting object bounding box
        NodeType first_geometry_low_node, first_geometry_high_node;
        r_edge_1.BoundingBox(first_geometry_low_node, first_geometry_high_node);

        // Creating OBB
        array_1d<array_1d<double, 3>, 2> first_direction_vector;
        noalias(first_direction_vector[0]) = first_geometry_high_node - first_geometry_low_node;
        const double norm_first_direction_vector = norm_2(first_direction_vector[0]);
        first_direction_vector[0] /= norm_first_direction_vector;
        first_direction_vector[1][0] = first_direction_vector[0][1];
        first_direction_vector[1][1] = - first_direction_vector[0][0];
        first_direction_vector[1][2] = 0.0;
        const array_1d<double, 3> first_center_point = 0.5 * (first_geometry_low_node.Coordinates() + first_geometry_high_node.Coordinates());
        array_1d<double, 2> first_half_distances;
        first_half_distances[0] = 0.5 * norm_first_direction_vector + mBoundingBoxFactor;
        first_half_distances[1] = mBoundingBoxFactor;
        OBB<2> first_obb(first_center_point, first_direction_vector, first_half_distances);

        // Second geometry
        for (std::size_t i_2 = 0; i_2 < number_of_edges_2; ++i_2) {
            auto& r_edge_2 = *(r_edges_2.begin() + i_2);

            // Check the intersection of each edge of the object bounding box against the intersecting object bounding box
            NodeType second_geometry_low_node, second_geometry_high_node;
            r_edge_2.BoundingBox(second_geometry_low_node, second_geometry_high_node);

            // Creating OBB
            array_1d<array_1d<double, 3>, 2> second_direction_vector;
            noalias(second_direction_vector[0]) = second_geometry_high_node - second_geometry_low_node;
            const double norm_second_direction_vector = norm_2(second_direction_vector[0]);
            second_direction_vector[0] /= norm_second_direction_vector;
            second_direction_vector[1][0] = second_direction_vector[0][1];
            second_direction_vector[1][1] = - second_direction_vector[0][0];
            second_direction_vector[1][2] = 0.0;
            const array_1d<double, 3> second_center_point = 0.5 * (second_geometry_low_node.Coordinates() + second_geometry_high_node.Coordinates());
            array_1d<double, 2> second_half_distances;
            second_half_distances[0] = 0.5 * norm_second_direction_vector + mBoundingBoxFactor;
            second_half_distances[1] = mBoundingBoxFactor;
            OBB<2> second_obb(second_center_point, second_direction_vector, second_half_distances);

            // Computing intersection OBB
            if (first_obb.HasIntersection(second_obb)){
                return true;
            }
        }
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::HasIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
//     // The local dimensions
//     const std::size_t local_dimension_1 = rFirstGeometry.LocalSpaceDimension();
//     const std::size_t local_dimension_2 = rSecondGeometry.LocalSpaceDimension();
//
//     // The faces
//     PointerVector<GeometryType> r_faces_1 = rFirstGeometry.Faces();
//     const std::size_t number_of_faces_1 = (local_dimension_1 < 2) ? 1 : r_faces_1.size();
//
//     PointerVector<GeometryType> r_faces_2 = rSecondGeometry.Faces();
//     const std::size_t number_of_faces_2 = (local_dimension_2 < 2) ? 1 : r_faces_2.size();
//
//     // First geometry
//     for (std::size_t i_1 = 0; i_1 < number_of_faces_1; ++i_1) {
//         auto& r_face_1 = *(r_faces_1.begin() + i_1);
//
//         // Check the intersection of each face of the object bounding box against the intersecting object bounding box
//         NodeType first_geometry_low_node, first_geometry_high_node;
//         r_face_1.BoundingBox(first_geometry_low_node, first_geometry_high_node);
//
//         // Creating OBB
//         array_1d<double, 3> first_direction_vector = first_geometry_high_node - first_geometry_low_node;
//         const double norm_first_direction_vector = norm_2(first_direction_vector);
//         first_direction_vector /= norm_first_direction_vector;
//         const array_1d<double, 3> first_center_point = 0.5 * (first_geometry_low_node.Coordinates() + first_geometry_high_node.Coordinates());
//         array_1d<double, 3> first_half_distances;
//         first_half_distances[0] = 0.5 * norm_first_direction_vector + mBoundingBoxFactor;
//         first_half_distances[1] = mBoundingBoxFactor;
//         first_half_distances[2] = mBoundingBoxFactor;
//         OBB<3> first_obb(first_center_point, first_direction_vector, first_half_distances);
//
//         // Second geometry
//         for (std::size_t i_2 = 0; i_2 < number_of_faces_2; ++i_2) {
//             auto& r_face_2 = *(r_faces_2.begin() + i_2);
//
//             // Check the intersection of each face of the object bounding box against the intersecting object bounding box
//             NodeType second_geometry_low_node, second_geometry_high_node;
//             r_face_2.BoundingBox(second_geometry_low_node, second_geometry_high_node);
//
//             // Creating OBB
//             array_1d<double, 3> second_direction_vector = second_geometry_high_node - second_geometry_low_node;
//             const double norm_second_direction_vector = norm_2(second_direction_vector);
//             second_direction_vector /= norm_second_direction_vector;
//             const array_1d<double, 3> second_center_point = 0.5 * (second_geometry_low_node.Coordinates() + second_geometry_high_node.Coordinates());
//             array_1d<double, 3> second_half_distances;
//             second_half_distances[0] = 0.5 * norm_second_direction_vector + mBoundingBoxFactor;
//             second_half_distances[1] = mBoundingBoxFactor;
//             second_half_distances[2] = mBoundingBoxFactor;
//             OBB<3> second_obb(second_center_point, second_direction_vector, second_half_distances);
//
//             // Computing intersection OBB
//             if (first_obb.HasIntersection(second_obb)){
//                 return true;
//             }
//         }
//     }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
Parameters FindIntersectedGeometricalObjectsWithOBBProcess<TEntity>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "first_model_part_name"  : "",
        "second_model_part_name" : "",
        "bounding_box_factor"    : -1.0
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class FindIntersectedGeometricalObjectsWithOBBProcess<Condition>;
template class FindIntersectedGeometricalObjectsWithOBBProcess<Element>;

}  // namespace Kratos.