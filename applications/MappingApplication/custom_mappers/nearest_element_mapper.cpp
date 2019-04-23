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
#include "nearest_element_mapper.h"
#include "mapping_application_variables.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos {

typedef std::size_t IndexType;
typedef std::size_t SizeType;

namespace {
bool ComputeProjection()

template<class TGeometryType>
void ComputeProjection(TGeometryType& rGeometry,
                       const Point& rPointToProject,
                       const double LocalCoordTol,
                       Vector& rShapeFunctionValues,
                       std::vector<int>& rEquationIds,
                       double& rDistance,
                       const ComputeApproximation=true)
{
    const SizeType num_nodes = p_geom->PointsNumber();
    const auto geom_family = p_geom->GetGeometryFamily();
    MapperUtilities::PairingIndex pairing_index = MapperUtilities::PairingIndex::Unspecified;
    std::vector<int> eq_ids;
    for (const auto& r_point : p_geom->Points()) {
        KRATOS_DEBUG_ERROR_IF_NOT(r_point.Has(INTERFACE_EQUATION_ID));
        eq_ids.push_back(r_point.GetValue(INTERFACE_EQUATION_ID));
    }

    if (geom_family == GeometryData::Kratos_Linear && num_nodes == 2) { // linear line
        pairing_index = MapperUtilities::ProjectOnLine(*p_geom, point_to_proj, LocalCoordTol, rShapeFunctionValues, eq_ids, proj_dist, ComputeApproximation);
    } else if ((geom_family == GeometryData::Kratos_Triangle      && num_nodes == 3) || // linear triangle
               (geom_family == GeometryData::Kratos_Quadrilateral && num_nodes == 4)) { // linear quad
        pairing_index = MapperUtilities::ProjectOnSurface(*p_geom, point_to_proj, LocalCoordTol, rShapeFunctionValues, eq_ids, proj_dist, ComputeApproximation);
    } else if (geom_family == GeometryData::Kratos_Tetrahedra ||
               geom_family == GeometryData::Kratos_Prism ||
               geom_family == GeometryData::Kratos_Hexahedra) { // Volume projection
        pairing_index = MapperUtilities::ProjectIntoVolume(*p_geom, point_to_proj, LocalCoordTol, rShapeFunctionValues, eq_ids, proj_dist, ComputeApproximation);
        // is_inside = MapperUtilities::ProjectIntoVolume(*p_geom, point_to_proj, local_coords, proj_dist);
    } else if (ComputeApproximation) {
        KRATOS_ERROR << "TODO Philipp implement a nearest neighbor here!" << std::endl;
    }

    const std::size_t num_values = shape_function_values.size();
    KRATOS_DEBUG_ERROR_IF_NOT(num_values == eq_ids.size()) << "Number of equation-ids is not the same as the number of ShapeFunction values!" << std::endl;

    return pairing_index;
}


}

void NearestElementInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject,
                                                      const double NeighborDistance)
{
    ComputeProjection();




    const auto p_geom = rInterfaceObject.pGetBaseGeometry();
    const SizeType num_nodes = p_geom->PointsNumber();

    const auto geom_family = p_geom->GetGeometryFamily();

    bool is_inside;
    double proj_dist;
    array_1d<double, 3> local_coords;

    const Point point_to_proj(this->Coordinates());

    // select projection depending on type of geometry
    if (geom_family == GeometryData::Kratos_Linear && num_nodes == 2) { // linear line
        Point projected_point;
        proj_dist = GeometricalProjectionUtilities::FastProjectOnLine(*p_geom, point_to_proj, projected_point);
        is_inside = p_geom->IsInside(projected_point, local_coords, 1e-14);
    } else if ((geom_family == GeometryData::Kratos_Triangle      && num_nodes == 3) || // linear triangle
               (geom_family == GeometryData::Kratos_Quadrilateral && num_nodes == 4)) { // linear quad
        Point projected_point;
        proj_dist = GeometricalProjectionUtilities::FastProjectOnGeometry(*p_geom, point_to_proj, projected_point);
        is_inside = p_geom->IsInside(projected_point, local_coords, 1e-14);
    } else if (geom_family == GeometryData::Kratos_Tetrahedra ||
               geom_family == GeometryData::Kratos_Prism ||
               geom_family == GeometryData::Kratos_Hexahedra) { // Volume projection
        is_inside = MapperUtilities::ProjectIntoVolume(*p_geom, point_to_proj, local_coords, proj_dist);
    } else {
        is_inside = false;
        KRATOS_WARNING_ONCE("NearestElementMapper") << "Unsupported type of geometry,"
            << "trying to use an approximation (Nearest Neighbor)" << std::endl;
    }

    Vector shape_function_values;

    // if it is closer, then we update the members to make this geometry the closest projection
    if (is_inside && proj_dist < mClosestProjectionDistance) {
        SetLocalSearchWasSuccessful();
        mClosestProjectionDistance = proj_dist;
        mShapeFunctionValues.clear();
        mNodeIds.clear();

        p_geom->ShapeFunctionsValues(shape_function_values, local_coords);
        KRATOS_DEBUG_ERROR_IF_NOT(shape_function_values.size() == num_nodes)
            << "Number of SFs is different from number of nodes!" << std::endl;

        if (mShapeFunctionValues.size() != num_nodes) mShapeFunctionValues.resize(num_nodes);
        if (mNodeIds.size() != num_nodes)             mNodeIds.resize(num_nodes);
        for (IndexType i=0; i<num_nodes; ++i) {
            mShapeFunctionValues[i] = shape_function_values[i];
            KRATOS_DEBUG_ERROR_IF_NOT((*p_geom)[i].Has(INTERFACE_EQUATION_ID))
                << "Node #" << (*p_geom)[i].Id() << " does not have an Interface Id!\n" << (*p_geom)<< "\n"
                <<  (*p_geom)[i] << std::endl;
            mNodeIds[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
        }
    }
}

void NearestElementInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject,
                                                                      const double NeighborDistance)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();
    const Point point_to_proj(this->Coordinates());
    Vector shape_function_values;
    double proj_dist;
    std::vector<int> eq_ids;
    for (const auto& r_point : p_geom->Points()) {
        KRATOS_DEBUG_ERROR_IF_NOT(r_point.Has(INTERFACE_EQUATION_ID));
        eq_ids.push_back(r_point.GetValue(INTERFACE_EQUATION_ID));
    }



    const std::size_t num_values = shape_function_values.size();
    KRATOS_DEBUG_ERROR_IF_NOT(num_values == eq_ids.size()) << "Number of equation-ids is not the same as the number of ShapeFunction values!" << std::endl;

    if (pairing_index > mPairingIndex || (pairing_index == mPairingIndex && proj_dist < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index;
        mClosestProjectionDistance = proj_dist;
        mNodeIds = eq_ids;

        if (mShapeFunctionValues.size() != num_values) mShapeFunctionValues.resize(num_values);
        for (std::size_t i=0; i<num_values; ++i) {
            mShapeFunctionValues[i] = shape_function_values[i];
        }
    }

    SetIsApproximation();
}

void NearestElementInterfaceInfo::SaveSearchResult()
{

}


void NearestElementLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() > 0) {
        double distance;
        double min_distance = std::numeric_limits<double>::max();
        int found_idx = -1;
        for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
            // the approximations will be processed in the next step if necessary
            if (!mInterfaceInfos[i]->GetIsApproximation()) {
                mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);
                if (distance < min_distance) {
                    min_distance = distance;
                    found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
                }
            }
        }

        if (found_idx == -1) { // no valid projection exists => using an approximation
            int int_pairing_index;
            MapperUtilities::PairingIndex pairing_index;
            for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
                // now the approximations are being checked
                if (mInterfaceInfos[i]->GetIsApproximation()) {
                    mInterfaceInfos[i]->GetValue(int_pairing_index, MapperInterfaceInfo::InfoType::Dummy);
                    pairing_index = (MapperUtilities::PairingIndex)int_pairing_index;
                    mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);

                    if (pairing_index > mPairingIndex || (pairing_index == mPairingIndex && distance < min_distance)) {
                        mPairingIndex = pairing_index;
                        min_distance = distance;
                        found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                        rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
                    }
                }
            }
        }

        KRATOS_ERROR_IF(found_idx == -1) << "Not even an approximation is found, this should not happen!"
            << std::endl; // TODO should thi sbe an error?

        std::vector<double> sf_values;

        mInterfaceInfos[found_idx]->GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
            rLocalMappingMatrix.resize(1, sf_values.size(), false);
        }
        for (IndexType i=0; i<sf_values.size(); ++i) {
            rLocalMappingMatrix(0,i) = sf_values[i];
        }

        mInterfaceInfos[found_idx]->GetValue(rOriginIds, MapperInterfaceInfo::InfoType::Dummy);

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string NearestElementLocalSystem::PairingInfo(const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestElementLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) {// TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
        if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
            mpNode->SetValue(PAIRING_STATUS, (int)mPairingIndex);
        }
    }
    return buffer.str();
}

}  // namespace Kratos.
