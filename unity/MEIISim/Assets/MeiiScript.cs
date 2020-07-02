using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

using fts;

public class MeiiScript : MonoBehaviour {

    [Header("Positions")]
    public float qe;
    public float qf;
    public float l1;
    public float l2;
    public float l3;
    public float theta1;
    public float theta2;
    public float theta3;
    public float px;
    public float py;
    public float pz;
    public float alpha;
    public float beta;
    public float gamma;

    [Header("Inertial Frame References")]
    public GameObject Elbow;
    public GameObject Forearm;
    public GameObject Rail1;
    public GameObject Rail2;
    public GameObject Rail3;
    public GameObject Slider1;
    public GameObject Slider2;
    public GameObject Slider3;
    public GameObject PlatformY;
    public GameObject PlatformZ;
    public GameObject PlatformX;

    private Vector3 Elbow_Zero;
    private Vector3 Forearm_Zero;
    private Vector3 Rail1_Zero;
    private Vector3 Rail2_Zero;
    private Vector3 Rail3_Zero;
    private Vector3 Slider1_Zero;
    private Vector3 Slider2_Zero;
    private Vector3 Slider3_Zero;
    private Vector3 PlatformY_Zero;
    private Vector3 PlatformZ_Zero;
    private Vector3 PlatformX_Zero;

    double[] qs = new double[14];

    // const string import_module = "virtual_meii";
    const string import_module = "meii_model";

    // bool printed = false;

    // Use this for initialization
    void Start () {
        // Elbow_Zero     = Elbow.transform.localEulerAngles;
        // Forearm_Zero   = Forearm.transform.localEulerAngles;
        // Rail1_Zero     = Rail1.transform.localPosition;
        // Rail2_Zero     = Rail2.transform.localPosition;
        // Rail3_Zero     = Rail3.transform.localPosition;
        // Slider1_Zero   = Slider1.transform.localEulerAngles;
        // Slider2_Zero   = Slider2.transform.localEulerAngles;
        // Slider3_Zero   = Slider3.transform.localEulerAngles;
        // PlatformY_Zero = PlatformY.transform.localEulerAngles;
        // PlatformZ_Zero = PlatformZ.transform.localEulerAngles;
        // PlatformX_Zero = PlatformX.transform.localEulerAngles;
        // Dll.stop();
        Dll.start();
	}

	// Update is called once per frame
	void Update () {
        Dll.get_positions(qs);
        qe     = (float)qs[0];
        qf     = (float)qs[1];
        l1     = (float)qs[2];
        l2     = (float)qs[3];
        l3     = (float)qs[4];
        theta1 = (float)qs[5];
        theta2 = (float)qs[6];
        theta3 = (float)qs[7];
        px     = (float)qs[8];
        py     = (float)qs[9];
        pz     = (float)qs[10];
        alpha  = (float)qs[12];
        beta   = (float)qs[13];
        gamma  = (float)qs[11];

        // foreach (var q in qs){
        //     Debug.Log(q);
        // }

        // Debug.Log(qe);
        // Debug.Log(qf);
        // Debug.Log(l1);
        // Debug.Log(l2);
        // Debug.Log(l3);
        // Debug.Log(theta1);
        // Debug.Log(theta2);
        // Debug.Log(theta3);
        // Debug.Log(px);
        // Debug.Log(py);
        // Debug.Log(pz);
        // Debug.Log(alpha);
        // Debug.Log(beta);
        // Debug.Log(gamma);
        
        Elbow.transform.localEulerAngles     = new Vector3(0, 0, -qe*Mathf.Rad2Deg);
        Forearm.transform.localEulerAngles   = new Vector3(qf*Mathf.Rad2Deg, 0, 0);
        Rail1.transform.localPosition        = new Vector3(-l1, Rail1_Zero.y, Rail1_Zero.z);
        Rail2.transform.localPosition        = new Vector3(-l2, Rail2_Zero.y, Rail2_Zero.z);
        Rail3.transform.localPosition        = new Vector3(-l3, Rail3_Zero.y, Rail3_Zero.z);
        Slider1.transform.localEulerAngles   = new Vector3( 5.415406f, 0, (Mathf.PI/2-theta1)*Mathf.Rad2Deg);
        Slider2.transform.localEulerAngles   = new Vector3(-114.5846f, 0, (Mathf.PI/2-theta2)*Mathf.Rad2Deg);
        Slider3.transform.localEulerAngles   = new Vector3( 125.4154f, 0, (Mathf.PI/2-theta3)*Mathf.Rad2Deg);
        PlatformY.transform.localPosition    = new Vector3(-px, py, pz);
        PlatformY.transform.localEulerAngles = new Vector3(PlatformY_Zero.x, -alpha*Mathf.Rad2Deg, PlatformY_Zero.z);
        PlatformZ.transform.localEulerAngles = new Vector3(PlatformZ_Zero.x, PlatformZ_Zero.y, -beta*Mathf.Rad2Deg);
        PlatformX.transform.localEulerAngles = new Vector3(gamma*Mathf.Rad2Deg, PlatformX_Zero.y, PlatformX_Zero.z);
    
        if (Input.GetKeyDown(KeyCode.R))
        {
            Dll.stop();
            Dll.start();
        }
    }

    void OnApplicationQuit() {
        Dll.stop();
    }

    // Dll Imports
    public class Dll {
        [DllImport(import_module)] 
        public static extern void start();
        [DllImport(import_module)] 
        public static extern void stop();
        // [DllImport(import_module)] 
        // public static extern void set_torques(double tau1, double tau2, double tau3, double tau4, double tau5);
        // [DllImport(import_module)]
        // public static extern void set_positions(double q1, double q2, double q3, double q4, double q5, double q6, double q7, double q8);
        // [DllImport(import_module)]
        // public static extern void set_velocities(double q1d, double q2d, double q3d, double q4d, double q5d, double q6d, double q7d, double q8d);
        [DllImport(import_module)]
        public static extern void get_positions(double[] positions);
    }

    // [PluginAttr("meii_model")]
    // public static class Dll
    // {        
    //     [PluginFunctionAttr("start")]
    //     public static Start start = null;
    //     public delegate void Start();
        
    //     [PluginFunctionAttr("stop")]
    //     public static Stop stop = null;
    //     public delegate void Stop();
        
    //     [PluginFunctionAttr("set_torques")]
    //     public static Set_torques set_torques = null;
    //     public delegate void Set_torques(double tau1, double tau2, double tau3, double tau4);
        
    //     [PluginFunctionAttr("set_positions")]
    //     public static Set_positions set_positions = null;
    //     public delegate void Set_positions(double q1, double q2, double q3, double q4, double q5, double q6, double q7);
        
    //     [PluginFunctionAttr("set_velocities")]
    //     public static Set_velocities set_velocities = null;
    //     public delegate void Set_velocities(double q1d, double q2d, double q3d, double q4d, double q5d, double q6d, double q7d);
        
    //     [PluginFunctionAttr("get_positions")]
    //     public static Get_positions get_positions = null;
    //     public delegate void Get_positions(double[] positions);
    // }
}

