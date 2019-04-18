//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "barycentric_mapper.h"
#include "mapping_application_variables.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{

namespace { // internal functions only used in this cpp

array_1d<double, 3> CreateArrayFromVector(const std::vector<double>& rVector,
                                          const std::size_t StartPosition) {
    array_1d<double,3> the_array;
    the_array[0] = rVector[StartPosition];
    the_array[1] = rVector[StartPosition+1];
    the_array[2] = rVector[StartPosition+2];
    return the_array;
}

void ComputeNeighborDistances(const array_1d<double,3>& rCoords,
                              const std::vector<double>& rNeighborCoords,
                              std::vector<double>& rDistances) {
    for (std::size_t i=0; i<rDistances.size(); ++i) {
        const auto neighbor_coords = CreateArrayFromVector(rNeighborCoords, i*3);
        rDistances[i] = MapperUtilities::ComputeDistance(rCoords, neighbor_coords);
    }
}

void InsertIfCloser(const array_1d<double,3>& rRefCoords,
                    const int CandidateEquationId,
                    const array_1d<double,3>& rCandidateCoords,
                    std::vector<int>& rNeighborIds,
                    std::vector<double>& rNeighborCoods) {
    const std::size_t num_interpolation_nodes = rNeighborIds.size();

    std::vector<double> neighbor_distances(num_interpolation_nodes, std::numeric_limits<double>::max());

    // compute the distances to the currently closest nodes based on the coordinates
    ComputeNeighborDistances(rRefCoords, rNeighborCoods, neighbor_distances);

    // compute distance to candidate
    const double candidate_distance = MapperUtilities::ComputeDistance(rRefCoords, rCandidateCoords);

    // check if the candidate is closer than the previously found nodes
    // save it if it is closer
    for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
        if (candidate_distance < neighbor_distances[i]) {
            rNeighborIds.insert(rNeighborIds.begin()+i, CandidateEquationId);
            rNeighborCoods.insert(rNeighborCoods.begin()+(i*3), std::begin(rCandidateCoords), std::end(rCandidateCoords));
            break;
        }
    }
    // resize is required because insert increases the size
    if (rNeighborIds.size() != num_interpolation_nodes) {
        rNeighborIds.resize(num_interpolation_nodes);
    }
    if (rNeighborCoods.size() != 3*num_interpolation_nodes) {
        rNeighborCoods.resize(3*num_interpolation_nodes);
    }
}

bool BarycentricInterpolateInEntity(const array_1d<double,3>& rRefCoords,
                                    const std::vector<double>& rCoordinates,
                                    Vector& rShapeFunctionValues) {

    bool is_inside;
    const std::size_t num_interpolation_nodes = rCoordinates.size()/3;

    array_1d<double, 3> local_coords;

    // TODO return false if points are coinciding and issue a warning, this should not happen!

    if (num_interpolation_nodes == 2) {
        Point::Pointer p_1(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 0)));
        Point::Pointer p_2(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 3)));
        Line2D2<Point> line(p_1, p_2);

        const Point point_to_proj(rRefCoords);
        double dummy = 0.0;

        is_inside = GeometricalProjectionUtilities::ProjectOnGeometry(line, point_to_proj, local_coords, dummy);
        if (is_inside) {
            line.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        }

    } else if (num_interpolation_nodes == 3) {
        Point::Pointer p_1(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 0)));
        Point::Pointer p_2(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 3)));
        Point::Pointer p_3(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 6)));
        Triangle3D3<Point> triangle(p_1, p_2, p_3);

        const Point point_to_proj(rRefCoords);
        double dummy = 0.0;

        // KRATOS_WATCH(rCoordinates)
        // KRATOS_WATCH(*p_1)
        // KRATOS_WATCH(*p_2)
        // KRATOS_WATCH(*p_3)

        // KRATOS_WATCH(triangle)

        is_inside = GeometricalProjectionUtilities::ProjectOnGeometry(triangle, point_to_proj, local_coords, dummy);
        if (is_inside) {
            triangle.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        }

    } else if (num_interpolation_nodes == 4) {
        Point::Pointer p_1(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 0)));
        Point::Pointer p_2(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 3)));
        Point::Pointer p_3(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 6)));
        Point::Pointer p_4(Kratos::make_shared<Point>(CreateArrayFromVector(rCoordinates, 9)));
        Tetrahedra3D4<Point> tetra(p_1, p_2, p_3, p_4);

        const Point point_to_proj(rRefCoords);
        double dummy = 0.0;

        is_inside = MapperUtilities::ProjectIntoVolume(tetra, point_to_proj, local_coords, dummy);
        if (is_inside) {
            tetra.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        }

    } else {
        KRATOS_ERROR << "Wrong number of interpolation nodes, this should not happen!" << std::endl;
    }

    KRATOS_WATCH(is_inside)

    return is_inside;
}

}

void BarycentricInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject,
                                                   const double NeighborDistance)
{
    SetLocalSearchWasSuccessful();

    const auto p_node = rInterfaceObject.pGetBaseNode();
    InsertIfCloser(
        Coordinates(),
        p_node->GetValue(INTERFACE_EQUATION_ID),
        p_node->Coordinates(),
        mNodeIds,
        mNeighborCoordinates);


    // KRATOS_WATCH(Coordinates())
    // KRATOS_WATCH(mNodeIds);
    // KRATOS_WATCH(mNeighborCoordinates);
    // std::cout << std::endl << std::endl;
}

void BarycentricLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() < 1) {
        ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
        return;
    }

    std::vector<int> node_ids;
    std::vector<double> neighbor_coods;

    // allocate final vectors, using the max possible size (in case of a volume interpolation)
    std::vector<int> final_node_ids(4);
    std::vector<double> final_neighbor_coords(12);
    std::fill(final_node_ids.begin(), final_node_ids.end(), -1);
    std::fill(final_neighbor_coords.begin(), final_neighbor_coords.end(), std::numeric_limits<double>::max());

    for (std::size_t i=0; i<mInterfaceInfos.size(); ++i) {
        mInterfaceInfos[i]->GetValue(node_ids, MapperInterfaceInfo::InfoType::Dummy);
        mInterfaceInfos[i]->GetValue(neighbor_coods, MapperInterfaceInfo::InfoType::Dummy);

        const std::size_t num_interpolation_nodes = node_ids.size();

        if (final_node_ids.size() != num_interpolation_nodes) {
            final_node_ids.resize(num_interpolation_nodes);
        }

        for (std::size_t j=0; j<num_interpolation_nodes; ++j) {
            InsertIfCloser(
                Coordinates(),
                node_ids[j],
                CreateArrayFromVector(neighbor_coods, j*3),
                final_node_ids,
                final_neighbor_coords);
        }
    }

    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
    rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

    if (std::find(final_node_ids.begin(), final_node_ids.end(), -1) == final_node_ids.end()) { // enough nodes found, trying to project/interpolate
        Vector shape_function_values;
        const bool is_inside = BarycentricInterpolateInEntity(Coordinates(), final_neighbor_coords, shape_function_values);
        if (is_inside) {
            rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
            if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != shape_function_values.size()) {
                rLocalMappingMatrix.resize(1, shape_function_values.size(), false);
            }
            for (std::size_t i=0; i<shape_function_values.size(); ++i) {
                rLocalMappingMatrix(0,i) = shape_function_values[i];
            }

            if (rOriginIds.size() != final_node_ids.size()) rOriginIds.resize(final_node_ids.size());
            for (std::size_t i=0; i<final_node_ids.size(); ++i) {
                rOriginIds[i] = final_node_ids[i];
            }

            return;
        }
    }

    // using the nearest neighbor as approximation
    rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
    if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) {
        rLocalMappingMatrix.resize(1, 1, false);
    }
    rLocalMappingMatrix(0,0) = 1.0;
    if (rOriginIds.size() != 1) rOriginIds.resize(1);
    rOriginIds[0] = final_node_ids[0];


    KRATOS_ERROR_IF(final_node_ids[0] == -1) << "Not even an approximation is found, this should not happen!" << std::endl; // TODO should this be an error?





    // KRATOS_WATCH(final_node_ids)
    // KRATOS_WATCH(final_neighbor_coords)


        // if enough_nodes_found:
        //     rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
        //     trying to project
        //     if this fails then also use nearest node
        // else:
        //     rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
        //     using nearest node
        //     if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) {
        //         rLocalMappingMatrix.resize(1, 1, false);
        //     }
        //     if (rOriginIds.size()      != 1) rOriginIds.resize(1);

        //     rLocalMappingMatrix(0,0) = 1.0;
        //     rOriginIds[0] = nearest_neighbor_id;

        // KRATOS_ERROR_IF(found_idx == -1) << "Not even an approximation is found, this should not happen!"
        //     << std::endl; // TODO should thi sbe an error?






        // //     // the approximations will be processed in the next step if necessary
        // //     if (!mInterfaceInfos[i]->GetIsApproximation()) {
        // //         mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);
        // //         if (distance < min_distance) {
        // //             min_distance = distance;
        // //             found_idx = static_cast<int>(i); // TODO explicit conversion needed?
        // //             rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
        // //         }
        // //     }
        // // }

        // // mInterfaceInfos[found_idx]->GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);
        // // reconstruct n closest nodes and create a geometry out of it
        // // then see if it is inside/outside


        // // rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

        // // if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) {
        // //     rLocalMappingMatrix.resize(1, 1, false);
        // // }
        // // if (rOriginIds.size()      != 1) rOriginIds.resize(1);
        // // if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

        // // int nearest_neighbor_id;
        // // double nearest_neighbor_distance;
        // // mInterfaceInfos[0]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
        // // mInterfaceInfos[0]->GetValue(nearest_neighbor_distance, MapperInterfaceInfo::InfoType::Dummy);

        // // for (std::size_t i=1; i<mInterfaceInfos.size(); ++i) {
        // //     // no check if this InterfaceInfo is an approximation is necessary
        // //     // bcs this does not exist for NearestNeighbor
        // //     double distance;
        // //     mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);

        // //     if (distance < nearest_neighbor_distance) {
        // //         nearest_neighbor_distance = distance;
        // //         mInterfaceInfos[i]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
        // //     }
        // // }

        // // rLocalMappingMatrix(0,0) = 1.0;
        // // rOriginIds[0] = nearest_neighbor_id;
        // // KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        // // rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

        // if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
        //     rLocalMappingMatrix.resize(1, sf_values.size(), false);
        // }
        // for (IndexType i=0; i<sf_values.size(); ++i) {
        //     rLocalMappingMatrix(0,i) = sf_values[i];
        // }

        // mInterfaceInfos[found_idx]->GetValue(rOriginIds, MapperInterfaceInfo::InfoType::Dummy);

}

std::string BarycentricLocalSystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "BarycentricLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    return buffer.str();
}

}  // namespace Kratos.
