bool test_boxshear_shearxy()
{
    printf("%s: %s: ",__FILE__,__FUNCTION__);
    ldf testers[]={0.0,0.0,0,0};
    ldf x[5]={-2., -1., 0., 1., 2.};
    ldf y[5]={0.0,0.0,0.0,0.0,0.0};
    ldf vx = .002;
    ldf dx[5]={vx,vx,vx,vx,vx};
    ldf dy[5]={0.0,0.0,0.0,0.0,0.0};
    md<2> sys(5);
    sys.set_rco(1.1);
    sys.set_ssz(1.1);
    sys.simbox.L[0]=5.0;
    sys.simbox.L[1]=5.0;
    // testing box matrix shear: shear boundaries perpendicular to x in y-direction
    // first index using ordinary PBC
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    vector<ldf> a={1.0,1.0};
    sys.add_typeinteraction(0,0,2,&a);
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    sys.index();
    sys.network.update=false;
    // now set shear of x boundary in y direction to -0.01
    sys.simbox.shear_boundary(1,0,-0.01);
    sys.integrator.method=INTEGRATOR::VVERLET;
    for (ui i = 0; i < 5; i++)
    {
        sys.particles[i].xp[0]=sys.particles[i].x[0]-sys.particles[i].dx[0]*sys.integrator.h;
        sys.particles[i].xp[1]=sys.particles[i].x[1]-sys.particles[i].dx[1]*sys.integrator.h;
    }
    for(ui h=0;h<400;h++)
    {
        for (ui i = 0; i < 5; i++) testers[0]+=sys.particles[i].x[0]; //fprintf(stdout,"%1.8Lf ",sys.particles[i].x[0]);
        for (ui i = 0; i < 5; i++) testers[1]+=sys.particles[i].x[0]; //fprintf(stdout,"%1.8Lf ",sys.particles[i].x[1]);
        for (ui i=0;i<2;i++) for (ui j=0;j<2;j++) testers[3]+=sys.simbox.Lshear[i][j]; //fprintf(stdout,"% 9.9Lf%c", sys.simbox.Lshear[i][j], j<1?' ':'\n');
        sys.timesteps(5000);
    }
    if(fabs(testers[0]+5.0)<=eps and fabs(testers[1]+5.0)<=eps and fabs(testers[2]-0.0)<=eps) return true;
    else return false;
}
