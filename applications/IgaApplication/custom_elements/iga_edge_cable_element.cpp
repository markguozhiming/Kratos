/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Anna Bauer
//                  Thomas Oberbichler
//                  Tobias Teschemacher
*/

// System includes
//#include "includes/define.h"
//#include "includes/variables.h"

// External includes

// Project includes
#include "iga_edge_cable_element.h"

namespace Kratos {              

Element::Pointer IgaEdgeCableElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const     
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaEdgeCableElement>(NewId, geometry,
        pProperties);
} 

void IgaEdgeCableElement::GetDofList(     
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void IgaEdgeCableElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetElementEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetElementEquationId(rResult, i, 2, DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void IgaEdgeCableElement::Initialize()
{
    mReferenceBaseVector = GetActualBaseVector();
}

//Definition von Base Vector 

IgaEdgeCableElement::Vector3 IgaEdgeCableElement::GetActualBaseVector() const
{
  // void CalculateTangent(
  //      const Matrix& rDN_De,
  //      const array_1d<double, 2>& Tangents);{
    
  //  return Tangents()
  //  }
    const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    const Vector& t = GetValue(TANGENTS);
    //const Vector& t = Tangent1;

    //const double t_u = t[0];
    //const double t_v = t[1];
    
    array_1d<double, 3> actual_base_vector = ZeroVector(3);

    IgaCurveOnSurfaceUtilities::CalculateTangent(
        GetGeometry(),
        DN_De,
        t,
        actual_base_vector
    );

KRATOS_WATCH(actual_base_vector)

    //for (std::size_t k = 0; k < NumberOfNodes(); k++) // k = Number of Nodes Cable
    //{
    //    actual_base_vector[0] += ( DN_De(k, 0) * t[0] + DN_De(k,1) * t[1] ) * GetGeometry()[k].X();        
    //    actual_base_vector[1] += ( DN_De(k, 0) * t[0] + DN_De(k,1) * t[1] ) * GetGeometry()[k].Y();        
    //    actual_base_vector[2] += ( DN_De(k, 0) * t[0] + DN_De(k,1) * t[1] ) * GetGeometry()[k].Z();
    //}

    return actual_base_vector;
}


void IgaEdgeCableElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

// get integration data

    const double& integration_weight = GetValue(INTEGRATION_WEIGHT);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    const Vector& t = GetValue(TANGENTS);

    // get properties

    const auto& properties = GetProperties();

    const double E = properties[YOUNG_MODULUS];
    const double A = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];

//KRATOS_WATCH(E)
//KRATOS_WATCH(A)
//KRATOS_WATCH(prestress)

    // compute base vectors

    const Vector3 actual_base_vector = GetActualBaseVector();

    const double reference_a = norm_2(mReferenceBaseVector);
    const double actual_a = norm_2(actual_base_vector);

    const double actual_aa = actual_a * actual_a;
    const double reference_aa = reference_a * reference_a;

    // green-lagrange strain

    const double e11_membrane = 0.0;//0.5 * (actual_aa - reference_aa);

    // normal force

    const double s11_membrane = prestress * A;//+ e11_membrane * A * E / reference_aa;

//KRATOS_WATCH(s11_membrane)

    //for (std::size_t k = 0; k < NumberOfDofs(); k++){
    //    const double dof_type_m[NumberOfDofs()];
    //    const std::size_t dof_type_m[ ] = { GetDofTypeIndex(k) };
    
    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        //const std::size_t dof_type_g[r] = {GetDofTypeIndex(r)};
        //const std::size_t dof_type_r = dof_type_m[r];
        const std::size_t dof_type_r = GetDofTypeIndex(r);
        const std::size_t shape_index_r = GetShapeIndex(r);
        //const double t_u = t[0];
        //const double t_v = t[1];
 
        const double epsilon_var_r = actual_base_vector[dof_type_r] *
            (shape_derivatives(shape_index_r, 0) * t[0] 
            + shape_derivatives(shape_index_r, 1) * t[1]) / reference_aa;
 
       if (ComputeLeftHandSide) {
            for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                const std::size_t dof_type_s = GetDofTypeIndex(s);
                const std::size_t shape_index_s = GetShapeIndex(s);
                const Vector& t = GetValue(TANGENTS);
                //const double t_u = t[0];
                //const double t_v = t[1];
 
                const double epsilon_var_s =
                    actual_base_vector[dof_type_s] *
                    (shape_derivatives(shape_index_s, 0) * t[0]
                    + shape_derivatives(shape_index_s, 1) * t[1])
                    / reference_aa;
 
                rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r *
                    epsilon_var_s;
 
                if (dof_type_r == dof_type_s) {
                    const double epsilon_var_rs =
                        (shape_derivatives(shape_index_r, 0) * t[0] + shape_derivatives(shape_index_r, 1) * t[1]) *
                        (shape_derivatives(shape_index_s, 0) * t[0] + shape_derivatives(shape_index_s, 1) * t[1]) / reference_aa;
 
                    rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs;
                }
            }
        }
        if (ComputeRightHandSide) {
            rRightHandSideVector[r] = -s11_membrane * epsilon_var_r;
        }
    }
    //}

//KRATOS_WATCH(reference_a)

    if (ComputeLeftHandSide) {
        rLeftHandSideMatrix *= reference_a * integration_weight;
    }

    if (ComputeRightHandSide) {
        rRightHandSideVector *= reference_a * integration_weight;
    }

//KRATOS_WATCH(rLeftHandSideMatrix)
//KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH("")
}

void IgaEdgeCableElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaEdgeCableElement\" #" << Id();
}

} // namespace Kratos