module imp_assembly_header
    use imp_assm_header
    implicit none
    type,public::sys_assembly!描述一个组件的特征和使用方法
      !private
      !real::fric  !摩擦因子
      type(hydraulic)::hydrau
      type(AssmGeom)::geom !assm_geom
      type(AssmMesh)::mesh !Assm_mesh
      type(material)::property !Assm_material 热物性和水力学参数
      type(th_boundary)::th_boundary !Assm_th_boundary
      type(AssmInit)::initdata
      type(confactor)::confactor_
      type(assmpow)::pow
      !real,allocatable::power(:,:) !Assm_power(zone,layer)
      !real,allocatable::fq_core(:,:)
      type(thermal)::Thermal  !pvt
    end type sys_assembly

end module imp_assembly_header
