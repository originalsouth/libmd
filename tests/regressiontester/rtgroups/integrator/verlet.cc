/* This test tests the Velocity Verlet integrator on correctness by integrating a particle from 0 to 1  */
bool test_integrator_verlet()
{
    ldf x[]={0.0};
    ldf y[]={0.0};
    ldf dx[]={0.1};
    ldf dy[]={0.1};
    md<2> test(1);
    test.network.update=false;
    test.integrator.h=1.0;
    test.integrator.method=INTEGRATOR::VVERLET;
    test.import_pos(x,y);
    test.import_vel(dx,dy);
    test.timesteps(10);
    test.export_pos(x,y);
    test.export_vel(dx,dy);
    if(fabs(dx[0]-0.1)<=eps and fabs(dy[0]-0.1)<=eps and fabs(x[0]-1.0)<=eps and fabs(y[0]-1.0)<=eps) test_success
    else test_fail
}
