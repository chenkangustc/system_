module imp_assembly_header
    use imp_assm_header
    implicit none
    type,public::sys_assembly!����һ�������������ʹ�÷���
      !private
      !real::fric  !Ħ������
      type(hydraulic)::hydrau
      type(AssmGeom)::geom !assm_geom
      type(AssmMesh)::mesh !Assm_mesh
      type(material)::property !Assm_material �����Ժ�ˮ��ѧ����
      type(th_boundary)::th_boundary !Assm_th_boundary
      type(AssmInit)::initdata
      type(confactor)::confactor_
      type(assmpow)::pow
      !real,allocatable::power(:,:) !Assm_power(zone,layer)
      !real,allocatable::fq_core(:,:)
      type(thermal)::Thermal  !pvt
    end type sys_assembly

end module imp_assembly_header
