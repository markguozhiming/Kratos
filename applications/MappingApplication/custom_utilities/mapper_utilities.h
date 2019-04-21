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

#if !defined(KRATOS_MAPPER_UTILITIES_H_INCLUDED)
#define  KRATOS_MAPPER_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
namespace MapperUtilities
{

enum class PairingIndex
{
    Volume_Inside   = -1,
    Volume_Outside  = -2,
    Surface_Inside  = -3,
    Surface_Outside = -4,
    Line_Inside     = -5,
    Line_Outside    = -6,
    Closest_Point   = -7,
    Unspecified     = -8
};

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Node<3> NodeType;

typedef Kratos::unique_ptr<MapperInterfaceInfo> MapperInterfaceInfoUniquePointerType;

typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
typedef std::vector<std::vector<MapperInterfaceInfoPointerType>> MapperInterfaceInfoPointerVectorType;

typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;
typedef Kratos::shared_ptr<MapperLocalSystemPointerVector> MapperLocalSystemPointerVectorPointer;


template< class TVarType >
static void FillFunction(const NodeType& rNode,
                         const TVarType& rVariable,
                         double& rValue)
{
    rValue = rNode.FastGetSolutionStepValue(rVariable);
}

template< class TVarType >
static void FillFunctionNonHist(const NodeType& rNode,
                                const TVarType& rVariable,
                                double& rValue)
{
    rValue = rNode.GetValue(rVariable);
}

template< class TVarType >
static std::function<void(const NodeType&, const TVarType&, double&)>
GetFillFunction(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::FROM_NON_HISTORICAL))
        return &FillFunctionNonHist<TVarType>;
    return &FillFunction<TVarType>;
}

template< class TVarType >
static void UpdateFunction(NodeType& rNode,
                           const TVarType& rVariable,
                           const double Value,
                           const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) = Value * Factor;
}

template< class TVarType >
static void UpdateFunctionWithAdd(NodeType& rNode,
                            const TVarType& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.FastGetSolutionStepValue(rVariable) += Value * Factor;
}

template< class TVarType >
static void UpdateFunctionNonHist(NodeType& rNode,
                            const TVarType& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) = Value * Factor;
}

template< class TVarType >
static void UpdateFunctionNonHistWithAdd(NodeType& rNode,
                            const TVarType& rVariable,
                            const double Value,
                            const double Factor)
{
    rNode.GetValue(rVariable) += Value * Factor;
}

template< class TVarType >
static std::function<void(NodeType&, const TVarType&, const double, const double)>
GetUpdateFunction(const Kratos::Flags& rMappingOptions)
{
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES) && rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNonHistWithAdd<TVarType>;
    if (rMappingOptions.Is(MapperFlags::ADD_VALUES))
        return &UpdateFunctionWithAdd<TVarType>;
    if (rMappingOptions.Is(MapperFlags::TO_NON_HISTORICAL))
        return &UpdateFunctionNonHist<TVarType>;
    return &UpdateFunction<TVarType>;
}

template< class TVectorType, class TVarType >
void UpdateSystemVectorFromModelPart(TVectorType& rVector,
                        ModelPart& rModelPart,
                        const TVarType& rVariable,
                        const Kratos::Flags& rMappingOptions)
{
    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        fill_fct(*(nodes_begin + i), rVariable, rVector[i]);
    }
}

template< class TVectorType, class TVarType >
void UpdateModelPartFromSystemVector(const TVectorType& rVector,
            ModelPart& rModelPart,
            const TVarType& rVariable,
            const Kratos::Flags& rMappingOptions)
{
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction<TVarType>(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        update_fct(*(nodes_begin + i), rVariable, rVector[i]);
    }
}

/**
* @brief Assigning INTERFACE_EQUATION_IDs to the nodes, with and without MPI
* This function assigns the INTERFACE_EQUATION_IDs to the nodes, which
* act as EquationIds for the MappingMatrix. This work with and without MPI,
* in MPI a ScanSum is performed with the local number of nodes
* @param rModelPartCommunicator The Modelpart-Communicator to be used
* @author Philipp Bucher
*/
void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator);

template<class TMapperLocalSystem>
void CreateMapperLocalSystemsFromNodes(const Communicator& rModelPartCommunicator,
                                       std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
{
    const std::size_t num_nodes = rModelPartCommunicator.LocalMesh().NumberOfNodes();
    const auto nodes_ptr_begin = rModelPartCommunicator.LocalMesh().Nodes().ptr_begin();

    if (rLocalSystems.size() != num_nodes) {
        rLocalSystems.resize(num_nodes);
    }

    #pragma omp parallel for
    for (int i = 0; i< static_cast<int>(num_nodes); ++i) {
        auto it_node = nodes_ptr_begin + i;
        rLocalSystems[i] = Kratos::make_unique<TMapperLocalSystem>((*it_node).get());
    }

    int num_local_systems = rLocalSystems.size(); // int bcs of MPI
    rModelPartCommunicator.SumAll(num_local_systems);

    KRATOS_ERROR_IF_NOT(num_local_systems > 0)
        << "No mapper local systems were created" << std::endl;
}

inline int ComputeNumberOfNodes(ModelPart& rModelPart)
{
    int num_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    rModelPart.GetCommunicator().SumAll(num_nodes); // Compute the sum among the partitions
    return num_nodes;
}

inline int ComputeNumberOfConditions(ModelPart& rModelPart)
{
    int num_conditions = rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
    rModelPart.GetCommunicator().SumAll(num_conditions); // Compute the sum among the partitions
    return num_conditions;
}

inline int ComputeNumberOfElements(ModelPart& rModelPart)
{
    int num_elements = rModelPart.GetCommunicator().LocalMesh().NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements); // Compute the sum among the partitions
    return num_elements;
}

template <class T1, class T2>
inline double ComputeDistance(const T1& rCoords1,
                              const T2& rCoords2)
{
    return std::sqrt( std::pow(rCoords1[0] - rCoords2[0] , 2) +
                      std::pow(rCoords1[1] - rCoords2[1] , 2) +
                      std::pow(rCoords1[2] - rCoords2[2] , 2) );
}

template<class TGeometryType>
PairingIndex ProjectOnLine(TGeometryType& rGeometry,
                           const Point& rPointToProject,
                           const double LocalCoordTol,
                           Vector& rShapeFunctionValues,
                           std::vector<int>& rEquationIds,
                           double& rDistance)
{
    Point projected_point;

    rDistance = GeometricalProjectionUtilities::FastProjectOnLine(rGeometry, rPointToProject, projected_point);

    array_1d<double, 3> local_coords;
    PairingIndex pairing_index;

    const bool is_inside = rGeometry.IsInside(projected_point, local_coords, 1e-14);

    if (is_inside) {
        pairing_index = PairingIndex::Line_Inside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
    } else if (!is_inside && std::abs(local_coords[0] < 1.0+LocalCoordTol)) {
        pairing_index = PairingIndex::Line_Outside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
    } else {
        // projection is ouside the line, searching the closest point
        pairing_index = PairingIndex::Closest_Point;
        const double dist_1 = ComputeDistance(rPointToProject, rGeometry[0]);
        const double dist_2 = ComputeDistance(rPointToProject, rGeometry[1]);

        if (dist_1 > dist_2) {
            rEquationIds[0] = rEquationIds[1];
            rDistance = dist_2;
        } else {
            rDistance = dist_1;
        }

        rEquationIds.resize(1);
        rShapeFunctionValues.resize(1);
        rShapeFunctionValues[0] = 1.0;
    }

    return pairing_index;
}

// TODO probably the eq-id handling can be improved after the interface-communication is changed
template<class TGeometryType>
PairingIndex ProjectOnSurface(TGeometryType& rGeometry,
                     const Point& rPointToProject,
                     const double LocalCoordTol,
                     Vector& rShapeFunctionValues,
                     std::vector<int>& rEquationIds,
                     double& rDistance)
{
    Point projected_point;

    rDistance = GeometricalProjectionUtilities::FastProjectOnGeometry(rGeometry, rPointToProject, projected_point);

    array_1d<double, 3> local_coords;
    PairingIndex pairing_index;

    const bool is_inside = rGeometry.IsInside(projected_point, local_coords, 1e-14);

    if (is_inside) {
        pairing_index = PairingIndex::Surface_Inside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
    } else if (!is_inside && std::abs(local_coords[0]) < 1.0+LocalCoordTol && std::abs(local_coords[1]) < 1.0+LocalCoordTol) {
        pairing_index = PairingIndex::Surface_Outside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
    } else { // inter-/extrapolation failed, trying to project on "subgeometries"
        pairing_index = PairingIndex::Unspecified;
        std::vector<int> final_eq_ids;

        auto edges = rGeometry.Edges();

        for (std::size_t i=0; i<edges.size(); ++i) {
            Vector edge_sf_values;
            std::vector<int> edge_eq_ids {rEquationIds[i], rEquationIds[(i+1)%edges.size()]};
            double edge_distance;

            const PairingIndex edge_index = ProjectOnLine(edges[i], rPointToProject, LocalCoordTol, edge_sf_values, edge_eq_ids, edge_distance);

            // check if the current edge gives a better result
            if (edge_index > pairing_index || (edge_index == pairing_index && edge_distance > rDistance)) {
                pairing_index = edge_index;
                rShapeFunctionValues = edge_sf_values;
                rDistance = edge_distance;
                final_eq_ids = edge_eq_ids;
            }
        }
        rEquationIds = final_eq_ids;
    }

    return pairing_index;
}


template<class TGeometryType>
bool ProjectIntoVolume(TGeometryType& rGeometry,
                       const Point& rPointToProject,
                       array_1d<double, 3>& rLocalCoords,
                       double& rDistance)
{
    bool is_inside = rGeometry.IsInside(rPointToProject, rLocalCoords);

    if (is_inside) {
        // Calculate Distance
        rDistance = ComputeDistance(rPointToProject, rGeometry.Center());
        rDistance /= rGeometry.Volume(); // Normalize Distance by Volume
    }

    return is_inside;
}

template <typename T>
inline double ComputeMaxEdgeLengthLocal(const T& rEntityContainer)
{
    double max_element_size = 0.0;
    // Loop through each edge of a geometrical entity ONCE
    for (const auto& r_entity : rEntityContainer) {
        for (std::size_t i = 0; i < (r_entity.GetGeometry().size() - 1); ++i) {
            for (std::size_t j = i + 1; j < r_entity.GetGeometry().size(); ++j) {
                double edge_length = ComputeDistance(r_entity.GetGeometry()[i].Coordinates(),
                                                        r_entity.GetGeometry()[j].Coordinates());
                max_element_size = std::max(max_element_size, edge_length);
            }
        }
    }
    return max_element_size;
}

inline double ComputeMaxEdgeLengthLocal(const ModelPart::NodesContainerType& rNodes)
{
    double max_element_size = 0.0;
    // TODO modify loop such that it loop only once over the nodes
    for (const auto& r_node_1 : rNodes) {
        for (const auto& r_node_2 : rNodes) {
            double edge_length = ComputeDistance(r_node_1.Coordinates(),
                                                    r_node_2.Coordinates());
            max_element_size = std::max(max_element_size, edge_length);
        }
    }
    return max_element_size;
}

double ComputeSearchRadius(ModelPart& rModelPart, int EchoLevel);

inline double ComputeSearchRadius(ModelPart& rModelPart1, ModelPart& rModelPart2, const int EchoLevel)
{
    double search_radius = std::max(ComputeSearchRadius(rModelPart1, EchoLevel),
                                    ComputeSearchRadius(rModelPart2, EchoLevel));

    KRATOS_INFO_IF("Mapper", EchoLevel > 0) << "Computed search-radius: "
        << search_radius << std::endl;

    return search_radius;
}

void CheckInterfaceModelParts(const int CommRank);

std::vector<double> ComputeLocalBoundingBox(ModelPart& rModelPart);

void ComputeBoundingBoxesWithTolerance(const std::vector<double>& rBoundingBoxes,
                                       const double Tolerance,
                                       std::vector<double>& rBoundingBoxesWithTolerance);

std::string BoundingBoxStringStream(const std::vector<double>& rBoundingBox);

bool PointIsInsideBoundingBox(const std::vector<double>& rBoundingBox,
                              const array_1d<double, 3>& rCoords);

void FillBufferBeforeLocalSearch(const MapperLocalSystemPointerVector& rMapperLocalSystems,
                                 const std::vector<double>& rBoundingBoxes,
                                 const SizeType BufferSizeEstimate,
                                 std::vector<std::vector<double>>& rSendBuffer,
                                 std::vector<int>& rSendSizes);

void CreateMapperInterfaceInfosFromBuffer(const std::vector<std::vector<double>>& rRecvBuffer,
                                          const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                          const int CommRank,
                                          MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer);

void FillBufferAfterLocalSearch(MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                const int CommRank,
                                std::vector<std::vector<char>>& rSendBuffer,
                                std::vector<int>& rSendSizes);

void AssignInterfaceInfosAfterRemoteSearch(const MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer,
                                           MapperLocalSystemPointerVectorPointer& rpMapperLocalSystems);

void DeserializeMapperInterfaceInfosFromBuffer(
    const std::vector<std::vector<char>>& rSendBuffer,
    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer);

/**
 * @class MapperInterfaceInfoSerializer
 * @ingroup MappingApplication
 * @brief Helper class to serialize/deserialize a vector containing MapperInterfaceInfos
 * @details This class serializes the vector containing the MapperInterfaceInfos (Shared Ptrs)
 * The goal of this class is to have a more efficient/faster implementation than the
 * one of the Serializer by avoiding the casting that is done in the serializer when pointers
 * are serialized
 * @TODO test the performance against the Serializer
 * @author Philipp Bucher
 */
class MapperInterfaceInfoSerializer
{
public:

    MapperInterfaceInfoSerializer(std::vector<MapperInterfaceInfoPointerType>& rMapperInterfaceInfosContainer,
                                  const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
        : mrInterfaceInfos(rMapperInterfaceInfosContainer)
        , mrpRefInterfaceInfo(rpRefInterfaceInfo->Create())
        { }

private:

    std::vector<MapperInterfaceInfoPointerType>& mrInterfaceInfos;
    MapperInterfaceInfoPointerType mrpRefInterfaceInfo;

    friend class Kratos::Serializer; // Adding "Kratos::" is nedded bcs of the "MapperUtilities"-namespace

    virtual void save(Kratos::Serializer& rSerializer) const;
    virtual void load(Kratos::Serializer& rSerializer);
};

}  // namespace MapperUtilities.

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_H_INCLUDED  defined
