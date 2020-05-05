﻿using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class MeiiScript : MonoBehaviour {

    [Header("Positions")]
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
    public GameObject Rail1;
    public GameObject Rail2;
    public GameObject Rail3;
    public GameObject Slider1;
    public GameObject Slider2;
    public GameObject Slider3;
    public GameObject PlatformY;
    public GameObject PlatformZ;
    public GameObject PlatformX;

    // [Header("Color Switching")]
    // public GameObject[] whites;
    // public GameObject[] colors;

    private Vector3 Rail1_Zero;
    private Vector3 Rail2_Zero;
    private Vector3 Rail3_Zero;
    private Vector3 Slider1_Zero;
    private Vector3 Slider2_Zero;
    private Vector3 Slider3_Zero;
    private Vector3 PlatformY_Zero;
    private Vector3 PlatformZ_Zero;
    private Vector3 PlatformX_Zero;

    double[] qs = new double[12];



    // Use this for initialization
    void Start () {
        Rail1_Zero     = Rail1.transform.localPosition;
        Rail2_Zero     = Rail2.transform.localPosition;
        Rail3_Zero     = Rail3.transform.localPosition;
        Slider1_Zero   = Slider1.transform.localEulerAngles;
        Slider2_Zero   = Slider2.transform.localEulerAngles;
        Slider3_Zero   = Slider3.transform.localEulerAngles;
        PlatformY_Zero = PlatformY.transform.localEulerAngles;
        PlatformZ_Zero = PlatformZ.transform.localEulerAngles;
        PlatformX_Zero = PlatformX.transform.localEulerAngles;
        Dll.start();
	}

	// Update is called once per frame
	void Update () {
        Dll.get_positions(qs);
        l1     = (float)qs[0];
        l2     = (float)qs[1];
        l3     = (float)qs[2];
        theta1 = (float)qs[3];
        theta2 = (float)qs[4];
        theta3 = (float)qs[5];
        px     = (float)qs[6];
        py     = (float)qs[7];
        pz     = (float)qs[8];
        alpha  = (float)qs[9];
        beta   = (float)qs[10];
        gamma  = (float)qs[11];
        
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

    /// Dll Imports
    public class Dll {
        [DllImport("meii_model")] 
        public static extern void start();
        [DllImport("meii_model")] 
        public static extern void stop();
        [DllImport("meii_model")] 
        public static extern void set_torques(double tau1, double tau2, double tau3);
        [DllImport("meii_model")]
        public static extern void set_positions(double q1, double q2, double q3, double q4, double q5, double q6);
        [DllImport("meii_model")]
        public static extern void set_velocities(double q1d, double q2d, double q3d, double q4d, double q5d, double q6d);
        [DllImport("meii_model")]
        public static extern void get_positions(double[] positions);
    }
}
