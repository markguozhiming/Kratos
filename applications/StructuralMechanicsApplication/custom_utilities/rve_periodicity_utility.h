//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED)
#define KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include <tuple>

// Project includes
#include "includes/define.h"
#include "includes/linear_master_slave_constraint.h"
#include "utilities/binbased_fast_point_locator_conditions.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class RVEPeriodicityUtility
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RVEPeriodicityUtility
    KRATOS_CLASS_POINTER_DEFINITION(RVEPeriodicityUtility);

    typedef std::tuple<std::vector<unsigned int>, std::vector<double>, Vector> DataTupletype;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RVEPeriodicityUtility(ModelPart &rDestinationModelPart) : mrModelPart(rDestinationModelPart) {}

    /// Destructor.
    virtual ~RVEPeriodicityUtility() {}

    ///@}
    ///@name Operators
    ///@{

    /** This function assign a pairing condition between two modelparts which contain two flat faces, parallel to 
     *  each other and separated by a distance rDistance.
     *  Note that this function should be called multiple times to pair the different faces in a box.
     * 
     *  @param rMasterModelPart master part to be paired
     *  @param rSlaveModelPart slave in the pairing
     *  @param rStrainTensor strain tensor which will be used in computing the pairing conditions
     *         the condition to be guaranteed will be that :    uslave = umaster + rStrainTensor * rDirection
     *  @param rDirection  a node with coordinates Xs on the slave, will be paired to the corresponding point with coordinates Xm on the master
     *         Xm will be computed as      Xm = Xs - rDirection
     */
    void AssignPeriodicity(ModelPart &rMasterModelPart,
                           ModelPart &rSlaveModelPart,
                           const Matrix &rStrainTensor,
                           const Vector &rDirection)
    {
        if (rMasterModelPart.Conditions().size() == 0)
            KRATOS_ERROR << "the master is expected to have conditions and it is empty" << std::endl;

        Vector translation = prod(rStrainTensor, rDirection);

        BinBasedFastPointLocatorConditions<3> bin_based_point_locator(rMasterModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        int max_search_results = 100;
        double search_tolerance = 1e-6;

        //construct auxiliary data structure to contain the master slave relation.
        //slave nodes must appear once, however a non-circular dependency is allowed between the masters
        for (unsigned int i = 0; i < rSlaveModelPart.Nodes().size(); ++i)
        {
            //search in which condition it falls
            auto i_node = rSlaveModelPart.NodesBegin() + i;

            Condition::Pointer p_host_cond;
            Vector N;
            array_1d<double, 3> transformed_slave_coordinates = i_node->Coordinates() - rDirection;

            // Finding the host element for this node
            const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(transformed_slave_coordinates, N, p_host_cond, max_search_results, search_tolerance);
            if (is_found)
            {
                const auto &geom = p_host_cond->GetGeometry();

                DataTupletype aux_data;

                auto &T = std::get<2>(aux_data);
                T = translation;

                auto &master_ids = std::get<0>(aux_data);
                auto &weights = std::get<1>(aux_data);
                for (unsigned int j = 0; j < geom.size(); ++j)
                {
                    master_ids.push_back(geom[j].Id());
                    weights.push_back(N[j]);
                }

                if (mAuxPairings.find(i_node->Id()) == mAuxPairings.end()) //this slave is not already present
                    mAuxPairings[i_node->Id()] = aux_data;
                else
                {
                    std::cout << "slave model part = " << rSlaveModelPart << std::endl;
                    std::cout << "master model part = " << rMasterModelPart << std::endl;
                    KRATOS_ERROR << "attempting to add twice the slave node with Id " << i_node->Id() << std::endl;
                }
            }
            else
            {
                KRATOS_ERROR << "counterpart not found for slave node " << i_node->Id() << std::endl;
            }
        }
    }

    /** this function finalizes the computation of the pairings. It can be called ONLY ONCE
     * @param rVariable is the value to which the pairing condition will be applied (needs to be a Variable with components)
     */
    void Finalize(const Variable<array_1d<double, 3>> &rVariable)
    {
        //get the components
        auto& rVar_x = KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(rVariable.Name() + "_X");
        auto& rVar_y = KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(rVariable.Name() + "_Y");
        auto& rVar_z = KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(rVariable.Name() + "_Z");

        for(auto& data : mAuxPairings)
        {
            auto &master_data = data.second;
            auto &master_ids = std::get<0>(master_data);
            auto &master_weights = std::get<1>(master_data);
            auto &T = std::get<2>(master_data);

            std::vector<unsigned int> final_master_ids;
            std::vector<double> final_master_weights;
            Vector final_T = T;

            for (unsigned int i = 0; i < master_ids.size(); ++i)
            {
                AppendIdsAndWeights(mAuxPairings, master_ids[i], master_weights[i], final_master_ids, final_master_weights, final_T);
            }

            //assign back the finalized pairings and weights to the data structure
            master_ids = final_master_ids;
            master_weights = final_master_weights;
            T = final_T;
        }

        //first assign master and slave all to false
        for (auto &node : mrModelPart.Nodes())
        {
            node.Set(SLAVE, false);
            node.Set(MASTER, false);
        }

        //compute the max id of the constraint
        unsigned int ConstraintId = 0;
        if (mrModelPart.MasterSlaveConstraints().size() != 0)
            ConstraintId = (mrModelPart.MasterSlaveConstraints().end() - 1)->Id();
        ConstraintId += 1;

        //define translation vector
        Vector xtranslation(1);
        Vector ytranslation(1);
        Vector ztranslation(1);

        for (const auto &data : mAuxPairings)
        {
            const unsigned int slave_id = data.first;
            const auto &master_data = data.second;
            auto &master_ids = std::get<0>(master_data);
            auto &master_weights = std::get<1>(master_data);
            auto &T = std::get<2>(master_data);

            //very useful for debugging. please do not remove
            // std::cout << slave_id << " - ";
            // for(auto& rMasterModelPart : master_ids)
            //     std::cout << rMasterModelPart << " ";
            // std::cout << " - ";
            // for(auto& w : master_weights)
            //     std::cout << w << " ";
            // std::cout << " - " << T << std::endl;

            //flag slave and master nodes
            mrModelPart.pGetNode(slave_id)->Set(SLAVE);
            for (auto &id : master_ids)
                mrModelPart.pGetNode(id)->Set(MASTER);

            

            //obtain the slave node
            auto pslave_node = mrModelPart.pGetNode(slave_id);

            //define relation matrix (same for the different components)
            Matrix relation_matrix(1, master_weights.size());
            for (unsigned int i = 0; i < relation_matrix.size2(); ++i)
                relation_matrix(0, i) = master_weights[i];

            xtranslation[0] = T[0];
            GenerateConstraint(ConstraintId,rVar_x, pslave_node, master_ids, relation_matrix, xtranslation);

            ytranslation[0] = T[1];
            GenerateConstraint(ConstraintId,rVar_y, pslave_node, master_ids, relation_matrix, ytranslation);

            ztranslation[0] = T[2];
            GenerateConstraint(ConstraintId,rVar_z, pslave_node, master_ids, relation_matrix, ztranslation);
        }
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "RVEPeriodicityUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "RVEPeriodicityUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{
    void AppendIdsAndWeights(
        std::map<unsigned int, DataTupletype> &rAux,
        const unsigned int MasterId,
        const double MasterWeight,
        std::vector<unsigned int> &rFinalMastersIds,
        std::vector<double> &rFinalMastersWeights,
        Vector &rFinalT)
    {
        if (std::abs(MasterWeight) > 1e-12) //discard nodes with negligible weight (note that weights sum to 1)
        {
            if (rAux.find(MasterId) == rAux.end()) //master is NOT also a slave
            {
                rFinalMastersIds.push_back(MasterId);
                rFinalMastersWeights.push_back(MasterWeight);
            }
            else //master also happens to be a slave
            {
                const auto &other_data = rAux[MasterId];
                const auto &other_master_ids = std::get<0>(other_data);
                const auto &other_master_weights = std::get<1>(other_data);
                const auto &other_T = std::get<2>(other_data);
                for (unsigned int j = 0; j < other_master_ids.size(); ++j)
                {
                    AppendIdsAndWeights(rAux, other_master_ids[j], MasterWeight * other_master_weights[j], rFinalMastersIds, rFinalMastersWeights, rFinalT);
                }

                rFinalT += MasterWeight * other_T;
            }
        }
    }

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart &mrModelPart;

    std::map<unsigned int, DataTupletype> mAuxPairings;

    ///@}
    ///@name Private Operators
    ///@{
    void GenerateConstraint(unsigned int& rConstraintId,
                            const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& rVar,
                            Node<3>::Pointer pSlaveNode,
                            const std::vector<unsigned int>& rMasterIds,
                            const Matrix& rRelationMatrix,
                            const Vector& rTranslationVector)
    {
        LinearMasterSlaveConstraint::DofPointerVectorType slave_dofs, master_dofs;
        slave_dofs.reserve(1);
        master_dofs.reserve(rMasterIds.size());

        slave_dofs.push_back(pSlaveNode->pGetDof(rVar));
        for (unsigned int i = 0; i < rMasterIds.size(); ++i)
            master_dofs.push_back(mrModelPart.pGetNode(rMasterIds[i])->pGetDof(rVar));
        mrModelPart.AddMasterSlaveConstraint(
            Kratos::make_shared<LinearMasterSlaveConstraint>(rConstraintId, master_dofs, slave_dofs, rRelationMatrix, rTranslationVector));
        rConstraintId += 1;
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RVEPeriodicityUtility &operator=(RVEPeriodicityUtility const &rOther) = delete;

    /// Copy constructor.
    RVEPeriodicityUtility(RVEPeriodicityUtility const &rOther) = delete;

    ///@}

}; // Class RVEPeriodicityUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                RVEPeriodicityUtility &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const RVEPeriodicityUtility &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED  defined