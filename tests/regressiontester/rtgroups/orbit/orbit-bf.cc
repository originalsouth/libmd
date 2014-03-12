/* This test tests if a two particle attrictive coulomb orbit is created properly in flat space */
bool test_orbit_orbit_bf()
{
    printf("%s: %s: ",__FILE__,__FUNCTION__);
    md<2> sys(2);
    ldf x[2]={-0.5,0.5};
    ldf y[2]={0.0,0.0};
    ldf dx[2]={0.0,0.0};
    ldf dy[2]={-0.5,0.5};
    ldf ans[4]={-0.486656996660,-0.080844601174,0.163903174469,-0.486480857881};
    sys.set_rco(10.0);
    sys.set_ssz(15.0);
    sys.simbox.L[0]=10.0;
    sys.simbox.L[1]=10.0;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.integrator.method=INTEGRATOR::VVERLET;
    sys.indexdata.method=INDEX::BRUTE_FORCE;
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    vector<ldf> a={-1.0};
    sys.add_typeinteraction(0,0,POT::POT_COULOMB,&a);
    sys.index();
    sys.network.update=false;
    for(ui h=0;h<500;h++) sys.timesteps(10);
    sys.export_pos(x,y);
    sys.export_vel(dx,dy);
    ldf sum_pos[2]={x[0]+x[1],y[0]+y[1]};
    ldf sum_vel[2]={dx[0]+dx[1],dy[0]+dy[1]};
    if(fabs(sum_pos[0])<=eps and fabs(sum_pos[1])<=eps and fabs(sum_vel[0])<=eps and fabs(sum_vel[1])<=eps and fabs(x[0]-ans[0])<=eps and fabs(y[0]-ans[1])<=eps and fabs(dx[0]-ans[2])<=eps and fabs(dy[0]-ans[3])<=eps) return true;
    else return false;
}
